"""Extract anchor and putative barcodes from FASTQ reads.

The extraction pipeline detects polyT tails and adapter placement, then scans for
barcode/UMI windows and classifies likely real anchors with a 2-component
negative-binomial mixture model.
"""

from Bio.Seq import reverse_complement
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter, defaultdict
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from contextlib import nullcontext
import logging
import os
from typing import Any

import edlib
from isal import igzip
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp
from scipy.stats import nbinom
from scipy.stats import norm
from sklearn.linear_model import LinearRegression

import config

logger = logging.getLogger(__name__)

# per-process worker state and helpers (used by ProcessPoolExecutor initializer)
_worker_processor = None


class _BarcodeBatchProcessor:
    """Worker-side barcode extraction processor.

    Instances are created once per process and reused for many record batches.
    """

    def __init__(self, polyt_min_length: int, adapter: str, barcode_len: int, umi_len: int, search_range: int,
                 min_qual: int, whitelist: set[str], write_upstream_region: bool = False):
        """
        Initialize the processor with parameters for barcode extraction.

        Args:
            polyt_min_length (int):       The minimum size required for a polyT(A) to be identified in the read
            adapter (str):                The adapter searched for in the read
            barcode_len (int):            The length of the cell barcode sequences
            umi_len (int):                The length of the UMI used
            search_range (int):           The range around the UMI length to look for exact matches
            min_qual (int):               The minimum average quality in a barcode for it to be considered
            whitelist (set):              The whitelist of cell barcodes
            write_upstream_region (bool): Whether to emit trimmed upstream sequence lines for output
        """
        self.barcode_len = barcode_len
        self.adapter = adapter
        self.umi_len = umi_len
        self.search_range = search_range
        self.min_qual = min_qual
        self.whitelist = whitelist
        self.write_upstream_region = write_upstream_region
        self.polyt_min_length = polyt_min_length
        self.buffer_dist = self.barcode_len + self.umi_len

    def process(self, records_batch: list[tuple]) -> dict:
        """
        Process a batch of FastQ records to extract barcodes.

        Args:
            records_batch (list): A list of tuples (title, sequence) for each read in the batch.

        Returns:
            dict:
                {
                    local_putative_lines: list[str],
                    local_upstream_region_lines: list[str],
                    local_putative_freq: Counter,
                    local_anchor_freq: Counter,
                    local_polya: int,
                    local_anchor: int,
                    local_quality_filtered_barcodes: int,
                    local_only_putative_reads: int,
                    local_processed: int,
                    local_both: int,
                    local_adapter: int
                }
        """
        local_putative_lines = []
        local_upstream_region_lines = []
        local_putative_freq = Counter()
        local_anchor_freq = Counter()
        local_polya = 0
        local_processed = 0
        local_anchor = 0
        local_quality_filtered_barcodes = 0
        local_only_putative_reads = 0
        local_both = 0
        local_adapter = 0

        for record in records_batch:
            matching_data = self._find_upstream(record[1], record[2])

            # counting
            if matching_data['f_polyt_score'] > 0 or matching_data['r_polyt_score'] > 0:
                local_polya += 1
            if matching_data['f_adapter_end'] or matching_data['r_adapter_end']:
                local_adapter += 1

            # check if both strands have a valid upstream, select the strand with anchors
            if matching_data['f_upstream_tail'] and matching_data['r_upstream_tail']:
                f_putative_barcodes, f_putative_umis, f_anchors, f_quality_filtered = self._extract_putatives(
                    matching_data['f_upstream_tail'], matching_data['f_upstream_qual'])
                r_putative_barcodes, r_putative_umis, r_anchors, r_quality_filtered = self._extract_putatives(
                    matching_data['r_upstream_tail'], matching_data['r_upstream_qual'])

                forward_len_diff = abs(len(matching_data['f_upstream_tail']) - (self.barcode_len + self.umi_len))
                reverse_len_diff = abs(len(matching_data['r_upstream_tail']) - (self.barcode_len + self.umi_len))

                # if only one orientation has an anchor, choose that
                if f_anchors and not r_anchors:
                    if f_putative_barcodes:
                        local_putative_lines.append(f'{record[0]}\t{",".join(f_putative_barcodes)}\t{",".join(f_putative_umis)}\t{matching_data["f_read"]}\t{matching_data["f_qual"]}\n')
                    for barcode in f_putative_barcodes:
                        local_putative_freq[barcode] += 1
                    for anchor in f_anchors:
                        local_anchor_freq[anchor] += 1
                        local_anchor += 1
                    local_quality_filtered_barcodes += f_quality_filtered

                elif not f_anchors and r_anchors:
                    if r_putative_barcodes:
                        local_putative_lines.append(f'{record[0]}\t{",".join(r_putative_barcodes)}\t{",".join(r_putative_umis)}\t{matching_data["r_read"]}\t{matching_data["r_qual"]}\n')
                    for barcode in r_putative_barcodes:
                        local_putative_freq[barcode] += 1
                    for anchor in r_anchors:
                        local_anchor_freq[anchor] += 1
                        local_anchor += 1
                    local_quality_filtered_barcodes += r_quality_filtered

                # choose the one with the proper length
                elif forward_len_diff < reverse_len_diff:
                    if f_putative_barcodes:
                        local_putative_lines.append(f'{record[0]}\t{",".join(f_putative_barcodes)}\t{",".join(f_putative_umis)}\t{matching_data["f_read"]}\t{matching_data["f_qual"]}\n')
                    for barcode in f_putative_barcodes:
                        local_putative_freq[barcode] += 1
                    for anchor in f_anchors:
                        local_anchor_freq[anchor] += 1
                        local_anchor += 1
                    local_quality_filtered_barcodes += f_quality_filtered
                elif reverse_len_diff < forward_len_diff:
                    if r_putative_barcodes:
                        local_putative_lines.append(f'{record[0]}\t{",".join(r_putative_barcodes)}\t{",".join(r_putative_umis)}\t{matching_data["r_read"]}\t{matching_data["r_qual"]}\n')
                    for barcode in r_putative_barcodes:
                        local_putative_freq[barcode] += 1
                    for anchor in r_anchors:
                        local_anchor_freq[anchor] += 1
                        local_anchor += 1
                    local_quality_filtered_barcodes += r_quality_filtered
                # else log that both strand could be valid
                else:
                    local_both += 1

            # check that only forward stand had a match
            elif matching_data['f_upstream_tail']:
                putative_barcodes, putative_umis, anchors, quality_filtered = self._extract_putatives(
                    matching_data['f_upstream_tail'], matching_data['f_upstream_qual'])

                if putative_barcodes:
                    local_putative_lines.append(f'{record[0]}\t{",".join(putative_barcodes)}\t{",".join(putative_umis)}\t{matching_data["f_read"]}\t{matching_data["f_qual"]}\n')
                for barcode in putative_barcodes:
                    local_putative_freq[barcode] += 1
                for anchor in anchors:
                    local_anchor_freq[anchor] += 1
                    local_anchor += 1
                if not anchors and putative_barcodes:
                    local_only_putative_reads += 1
                local_quality_filtered_barcodes += quality_filtered

            # check that only reverse strand had a match
            elif matching_data['r_upstream_tail']:
                putative_barcodes, putative_umis, anchors, quality_filtered = self._extract_putatives(
                    matching_data['r_upstream_tail'], matching_data['r_upstream_qual'])

                if putative_barcodes:
                    local_putative_lines.append(f'{record[0]}\t{",".join(putative_barcodes)}\t{",".join(putative_umis)}\t{matching_data["r_read"]}\t{matching_data["r_qual"]}\n')
                for barcode in putative_barcodes:
                    local_putative_freq[barcode] += 1
                for anchor in anchors:
                    local_anchor_freq[anchor] += 1
                    local_anchor += 1 
                if not anchors and putative_barcodes:
                    local_only_putative_reads += 1
                local_quality_filtered_barcodes += quality_filtered
            local_processed += 1

        return {
            'local_putative_lines': local_putative_lines,
            'local_upstream_region_lines': local_upstream_region_lines,
            'local_putative_freq': local_putative_freq,
            'local_anchor_freq': local_anchor_freq,
            'local_polya': local_polya,
            'local_anchor': local_anchor,
            'local_quality_filtered_barcodes': local_quality_filtered_barcodes,
            'local_only_putative_reads': local_only_putative_reads,
            'local_processed': local_processed,
            'local_both': local_both,
            'local_adapter': local_adapter
        }

    @staticmethod
    def _avg_qual(qual_str: str) -> float:
        """
        Determine the average quality of a Phred encoded quality string

        Args:
            qual_str (str): Phred encoded quality string

        Returns:
            float: average quality of the string
        """
        total = 0
        for qchar in qual_str:
            total += ord(qchar) - 33
        return total / len(qual_str)

    def _extract_putatives(self, upstream_tail: str, upstream_qual: str) -> tuple[list[str], list[str], list[str], int]:
        """Collect all putative barcode/UMI pairs and mark which are anchors.

        Args:
            upstream_tail (str): The sequence segment upstream of the polyT to search.
            upstream_qual (str): Quality string for upstream_tail.

        Returns:
            list[str]: putative barcode sequences
            list[str]: putative UMI sequences
            list[str]: putative anchor sequences
            int:       number of barcodes filtered for quality
            
        """
        putative_barcodes = []
        putative_umis = []
        anchors = []
        quality_filtered = 0

        expected_pos = 0
        upper_search = expected_pos + self.search_range

        for i in range(0, upper_search):
            checker = upstream_tail[i:i + self.barcode_len]
            checker_qual = upstream_qual[i:i + self.barcode_len]
            umi = upstream_tail[i + self.barcode_len:i + self.barcode_len + self.umi_len]
            if len(checker) != self.barcode_len:
                quality_filtered += 1
                continue
            if checker in self.whitelist and self._avg_qual(checker_qual) >= self.min_qual:
                anchors.append(checker)
            putative_barcodes.append(checker)
            # add T's to the UMI when too short from running the window upstream and including the identified polyT in the UMI
            while len(umi) != self.umi_len:
                umi = umi + 'T'
            putative_umis.append(umi)

        return putative_barcodes, putative_umis, anchors, quality_filtered

    def _adapter_search(self, seq: str) -> int | None:
        """
        Identify adapter sequence in a region upstream of polyT.

        Args:
            seq (str): upstream region to search

        Returns:
            int | None : End position of the adapter, or None if not found
        """
        max_mismatches = 2
        upstream_len = len(seq)

        aln = edlib.align(self.adapter, seq, mode='HW', task='locations', k=max_mismatches)
        edit_dist = aln['editDistance']

        # check if no alignment
        if edit_dist == -1:
            return None

        # select the alignment that yields the most accurate barcode region length
        end = None
        current_barcode_region_len = None
        for location in aln['locations']:
            if end == None:
                end = location[1]
                current_barcode_region_len = upstream_len - end
            elif abs(current_barcode_region_len  - self.buffer_dist) > abs((upstream_len - location[1]) - self.buffer_dist):
                end = location[1]
                current_barcode_region_len = upstream_len - end

        return end + 1

    def _polyt_search(self, seq: str) -> dict:
        """
        Find the location of the polyT in a read. Returns the position and length
        if it mets the minimum length

        Args:
            seq (str): the read

        Returns:
            dict:
                {
                    best_start: int, 
                    best_len: int, 
                    best_score: int
                }
        """
        # initialize
        end_trim = 2
        score = 0.0
        best_score = 0.0
        best_start = 0
        best_len = 0
        trailing_t = 0
        best_valid_end = -1  # last position where trailing T is strong
        mismatches = 0
        i = 0

        while i < len(seq):
            if score == 0:
                if seq[i:i + self.polyt_min_length] != "T" * self.polyt_min_length:
                    i += 1
                    continue

                current_start = i
                current_len = self.polyt_min_length
                score = self.polyt_min_length

                trailing_t = self.polyt_min_length
                best_valid_end = i + self.polyt_min_length - 1

                if score > best_score:
                    best_score = score
                    best_start = current_start
                    best_len = current_len

                i += self.polyt_min_length
            else:
                if seq[i] == "T":
                    score += 1
                    current_len += 1
                    trailing_t += 1
                else:
                    score -= 2
                    current_len += 1
                    trailing_t = 0
                    mismatches += 1

                # record valid endpoints when trailing T is strong enough
                if trailing_t >= end_trim:
                    best_valid_end = i

                if score <= 0:
                    score = 0
                    current_len = 0
                    trailing_t = 0
                    best_valid_end = -1
                    mismatches = 0
                else:
                    if best_valid_end != -1:
                        curr_len = best_valid_end - current_start + 1
                        if score > best_score:
                            best_score = score
                            best_start = current_start
                            best_len = curr_len

                i += 1

        return {
                'start': best_start, 
                'len': best_len, 
                'score': best_score
                }

    def _find_upstream(self, sequence: str, qual: str) -> dict:
        """
        Find and return the subsequence upstream of the polyT.

        Args:
            sequence (str): Full read sequence.
            qual (str):     Quality string for the read

        Returns:
            dict:
                {
                    f_polyt_start: int | None,
                    f_polyt_score: int,
                    r_polyt_start: int | None,
                    r_polyt_score: int,
                    f_adapter_end: int | None,
                    r_adapter_end: int | None,
                    f_upstream_tail: str | None,
                    f_upstream_qual: str | None,
                    r_upstream_tail: str | None,
                    r_upstream_qual: str | None,
                    f_read: str | None,
                    f_qual: str | None,
                    r_read: str | None,
                    r_qual: str | None
                }
        """
        # initialize
        reverse_seq = reverse_complement(sequence)
        reverse_qual = qual[::-1]
        buffer_bases = self.barcode_len + self.umi_len
        f_adapter_end, r_adapter_end = None, None
        returning = {'f_polyt_start': None,
                     'f_polyt_score': 0,
                     'r_polyt_start': None,
                     'r_polyt_score': 0,
                     'f_adapter_end': None,
                     'r_adapter_end': None,
                     'f_upstream_tail': None,
                     'f_upstream_qual': None,
                     'r_upstream_tail': None,
                     'r_upstream_qual': None,
                     'f_read': None,
                     'f_qual': None,
                     'r_read': None,
                     'r_qual': None
                     }

        # forward strand
        f_polyt = self._polyt_search(sequence[buffer_bases:250])
        f_polyt['start'] += buffer_bases
        if f_polyt['score'] > 0:
            f_adapter_end = self._adapter_search(sequence[:f_polyt['start']])
            returning['f_polyt_start'] = f_polyt['start']
            returning['f_polyt_score'] = f_polyt['score']
            returning['f_adapter_end'] = f_adapter_end

        # reverse strand
        r_polyt = self._polyt_search(reverse_seq[buffer_bases:250])
        r_polyt['start'] += buffer_bases
        if r_polyt['score'] > 0:
            r_adapter_end = self._adapter_search(reverse_seq[:r_polyt['start']])
            returning['r_polyt_start'] = r_polyt['start']
            returning['r_polyt_score'] = r_polyt['score']
            returning['r_adapter_end'] = r_adapter_end

        # check if strands hit
        if f_polyt['score'] > 0 and f_adapter_end:
            returning['f_upstream_tail'] = sequence[f_adapter_end:f_polyt['start'] + self.search_range]
            returning['f_upstream_qual'] = qual[f_adapter_end:f_polyt['start'] + self.search_range]
            returning['f_read'] = sequence[f_polyt['start'] + f_polyt['len']:]
            returning['f_qual'] = qual[f_polyt['start'] + f_polyt['len']:]

        if r_polyt['score'] > 0 and r_adapter_end:
            returning['r_upstream_tail'] = reverse_seq[r_adapter_end: r_polyt['start'] + self.search_range]
            returning['r_upstream_qual'] = reverse_qual[r_adapter_end: r_polyt['start'] + self.search_range]
            returning['r_read'] = reverse_seq[r_polyt['start'] + r_polyt['len']:]
            returning['r_qual'] = reverse_qual[r_polyt['start'] + r_polyt['len']:]

        return returning    


def _worker_init(polyt_min_length: int, adapter: str, barcode_len: int, umi_len: int, search_range: int, min_qual: int, whitelist: set[str],
                 write_upstream_region: bool) -> None:
    """
    Initializer executed in each worker process to set up a local processor.

    Args:
        polyt_min_length (int): The minimum size required for a polyT(A) to be identified in the read
        adatper (str):          The adapter to be searched for in the read
        barcode_len (int):      The length of the cell barcode sequences
        umi_len (int):          The length of the UMI used
        search_range (int):     The range around the UMI length to look for exact matches
        min_qual (int):         The minimum average quality in a barcode for it to be considered
        whitelist (set):        The whitelist of cell barcodes
        write_upstream_region (bool): Whether to emit trimmed upstream sequence lines for output
    """
    global _worker_processor
    _worker_processor = _BarcodeBatchProcessor(polyt_min_length, adapter, barcode_len, umi_len, search_range,
                                              min_qual, whitelist,
                                              write_upstream_region=write_upstream_region)


def _worker_process_batch(records_batch: list[tuple[str, str, str]]) -> dict[str, Any]:
    """
    Module-level entry point used by ProcessPoolExecutor workers.

    Args:
        records_batch (list): A list of tuples (title, sequence, qual) for each read in the batch.

    Returns:
        dict:
            {
                local_putative_lines: list[str],
                local_upstream_region_lines: list[str],
                local_putative_freq: Counter,
                local_anchor_freq: Counter,
                local_polya: int,
                local_anchor: int,
                local_quality_filtered_barcodes: int,
                local_only_putative_reads: int,
                local_processed: int,
                local_both: int,
                local_adapter: int
            }
    Raises:
        RuntimeError
    """
    if _worker_processor is None:
        raise RuntimeError(
            'Worker not initialized. This function must be run inside a worker process with _worker_init called.')
    return _worker_processor.process(records_batch)


def _safe_future_result(future: Future[dict[str, Any]]) -> dict[str, Any]:
    """
    Return future.result() ensuring the exception is logged with traceback.
    Worker exceptions are logged here then re-raised so callers still fail fast.

    Returns:
        dict:
            {
                local_putative_lines: list[str],
                local_upstream_region_lines: list[str],
                local_putative_freq: Counter,
                local_anchor_freq: Counter,
                local_polya: int,
                local_anchor: int,
                local_quality_filtered_barcodes: int,
                local_only_putative_reads: int,
                local_processed: int,
                local_both: int,
                local_adapter: int
            }
    Raises:
        Exception: From the worker
    """
    try:
        return future.result()
    except Exception:
        logger.exception('Worker task raised an exception (see traceback) — re-raising')
        raise


def _estimate_nb_params(vals: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """
    weighted method-of-moments estimate for negative binomial parameters

    Args:
        vals (np.array):    count data
        weights (np.array): weights assigned to data

    Returns:
        int: r value
        int: p value
    """
    weight_sum = weights.sum()
    if weight_sum <= 0:
        w = np.full(len(vals), 1.0 / len(vals), dtype=float)
    else:
        w = weights / weight_sum
    mean = np.sum(w * vals)
    var = np.sum(w * (vals - mean) ** 2)
    var = max(var, mean + 1e-3)  # don't allow it to be too low

    r = mean ** 2 / (var - mean)
    r = max(r, 1e-3)

    p = r / (r + mean)
    p = float(np.clip(p, 1e-6, 1 - 1e-6))

    return r, p


def _fit_nb_mixture_2component(anchor_dist: Counter[str], probability_threshold: float, output_csv: str | None = None,
                              expected_cells: int | None = None, prior_strength: int = 200, max_iter: int = 200,
                              tol: float = 1e-6) -> dict[str, Any]:
    """
    Deterministic 2-component NB mixture (noise + real).

    Args:
        anchor_counter (Counter):      Dict of anchor: frequency
        probability_threshold (float): Cutoff to consider a real from the real component
        expected_cells (int | None):   Number of cells expected to have been sequenced
        prior_strength (int):          Strength of belief of prior
        max_iter (int):                Number of EM iterations
        tol (float):                   Tolerance of the model
        output_csv (str):              Name of file to write probabilities to

    Returns:
        dict of structure:
            {
            probs_real: list[float],
            probs_noise: list[float],
            assignments: list[int],
            threshold: int,
            params: {
                noise: {
                    r: float,
                    p: float
                },
                real: {
                    r: float,
                    p: float
                },
                mixing_fraction: float,
            }
    """
    x = np.array(sorted(list(set(anchor_dist.values())), reverse=True), dtype=float)
    n = len(x)

    ranks = np.arange(n)

    real_cut = max(1, n // 2)

    z_noise = (ranks >= real_cut).astype(float)
    z_real = (ranks < real_cut).astype(float)

    if expected_cells is not None:
        expected_fraction = min(max(expected_cells / n, 1e-6), 1 - 1e-6)
        alpha = np.array([1.0 - expected_fraction, expected_fraction], dtype=float) * prior_strength
        pi = alpha / alpha.sum()
    else:
        alpha = np.ones(2, dtype=float)
        pi = np.array([z_noise.mean(), z_real.mean() ], dtype=float)
        pi = np.clip(pi, 1e-6, None)
        pi /= pi.sum()

    r_noise, p_noise = _estimate_nb_params(x, z_noise + 1e-3)
    r_real, p_real = _estimate_nb_params(x, z_real + 1e-3)

    prev_ll = -np.inf

    for _ in range(max_iter):

        log_probs = np.vstack([
            np.log(pi[0]) + nbinom.logpmf(x, r_noise, p_noise),
            np.log(pi[1]) + nbinom.logpmf(x, r_real, p_real),
        ])

        denom = logsumexp(log_probs, axis=0)
        gamma = np.exp(log_probs - denom)

        if expected_cells is not None:
            pi = (gamma.sum(axis=1) + alpha - 1) / (len(x) + alpha.sum() - len(alpha))
        else:
            pi = gamma.mean(axis=1)

        pi = np.clip(pi, 1e-6, None)
        pi /= pi.sum()

        r_noise, p_noise = _estimate_nb_params(x, gamma[0])
        r_real, p_real = _estimate_nb_params(x, gamma[1])

        # enforce ordering (real must have larger mean)
        mean_noise = r_noise * (1 - p_noise) / p_noise
        mean_real = r_real * (1 - p_real) / p_real

        if mean_noise > mean_real:
            r_noise, r_real = r_real, r_noise
            p_noise, p_real = p_real, p_noise
            gamma = gamma[::-1]
            pi = pi[::-1]

        ll = denom.sum()

        if np.abs(ll - prev_ll) < tol:
            break

        prev_ll = ll

    probs_real = gamma[1]

    assignments = (probs_real > probability_threshold).astype(int)

    real_counts = np.sort(x[assignments == 1])
    threshold = real_counts[0] if len(real_counts) > 0 else None

    if output_csv is not None:
        with open(output_csv, "w") as f:
            f.write("barcode,count,prob_noise,prob_real\n")

            sorted_items = sorted(anchor_dist.items(), key=lambda kv: kv[1], reverse=True)

            for (barcode, count), noise_prob, real_prob in zip(sorted_items, gamma[0], gamma[1]):
                f.write(f"{barcode},{count},{noise_prob:.6f},{real_prob:.6f}\n")

    return {
        "probs_real": probs_real,
        "probs_damage": np.zeros_like(probs_real),
        "probs_noise": gamma[0],
        "assignments": assignments,
        "threshold": threshold,
        "params": {
            "noise": {
                "r": r_noise,
                "p": p_noise
            },
            "real": {
                "r": r_real,
                "p": p_real
            },
            "mixing_fraction": float(pi[1]),
        },
    }


def _plot_anchor_probs(anchor_counter: Counter[str], probs_real: np.ndarray, probs_noise: np.ndarray, cutoff: int | None, probability_threshold: float,
                      outfile: str) -> None:
    """
    Plot:
    1) ranked anchor counts
    2) count vs probability

    Args:
        anchor_counter (Counter):      barcode -> count
        probs_real (array):            P(real | x) aligned to sorted counts (desc)
        cutoff (int):                  NB-derived cutoff
        probability_threshold (float): Probability cutoff for determining cells
        outfile (str):                 The name of the file to write probability information to
    """
    counts = np.array(sorted(list(set(anchor_counter.values())), reverse=True), dtype=float)

    real_threshold = np.sum(counts > cutoff)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.scatter(np.arange(len(counts)), counts, s=5)
    ax1.axvline(x=real_threshold, linestyle='--', color='red', label='NB cutoff')
    ax1.set_title('Anchor Frequency (Ranked)')
    ax1.set_xlabel('Ranked Index')
    ax1.set_ylabel('Frequency')
    ax1.legend()

    
    ax2.scatter(counts, probs_real, s=5, label='P(real)', color='tab:blue')
    ax2.scatter(counts, probs_noise, s=5, label='P(noise)', color='tab:green')
    ax2.set_title('P(real | count)')
    ax2.set_xlabel('Anchor Frequency')
    ax2.set_ylabel('Probability')
    ax2.axhline(y=probability_threshold, linestyle='--', color='gray', label='Probability Threshold')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def extract_barcodes(fastq: str, polyt_min_length: int, adapter: str, barcode_len: int, umi_len: int,
                     search_range: int, min_qual: int, expected_cells: int,
                     whitelist: set[str] | None = None, processes: int | None = None, batch_size: int = 10000,
                     upstream_region_tsv: str | None = None, probability_threshold: float = 0.85, max_iter: int = 200) -> tuple[dict[str, int], dict[str, int]]:
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):                     The path to the fastq file
        polyt_min_length (int):          The minimum size required for a polyT(A) to be identified in the read
        adapter (str):                   The adapter to be searched for in the read
        barcode_len (int):               The length of the cell barcode sequenceS
        umi_len (int):                   The length of the UMI used
        search_range (int):              The range around the UMI length to look for exact matches
        min_qual (int):                  The minimum average quality in a barcode for it to be considered
        expected_cells (int):            The number of cells that were expected to be sequenced
        whitelist (set):                 The whitelist of cell barcodes that are expected as a set for quick accession
        processes (int):                 Number of worker processes (defaults to cpu_count() - 1)
        batch_size (int):                Number of reads per worker task (default 1000)
        upstream_region_tsv (str):       Output path for optional TSV of search-window-to-polyT sequences
        probability_threshold (float):   The probability required from the NB mixture to be considered a true anchor.
        max_iter (int):                  The maximum number of iteration for EM.

    Returns:
        defaultdict: anchor barcode: frequency
        defaultdict: putative barcode: frequency
    """
    logger.section('BARCODE EXTRACTATION')
    logger.info('Beginning extraction')

    # initialize values
    putative_anchors_freq = Counter()
    putative_barcodes_freq = Counter()
    polya_found_reads = 0
    reads_with_putative_anchor_count = 0
    processed = 0
    total_quality_filtered_barcodes = 0
    only_putative_barcodes = 0
    both_strands_valid = 0
    adapter_found = 0
    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)
    
    # initialize file handles
    putative_barcode_file = os.path.join(config.GLOBAL_OUTDIR, '.putative_barcodes')
    anchor_freqs_file = os.path.join(config.GLOBAL_OUTDIR, '.anchor_freqs')
    probability_figure = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_anchor_frequency_thresholds_and_probability.png')
    probabilities_csv = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_anchor_probabilities.csv')

    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)

    # stream through fastq file and submit batches to worker processes (decompression/streaming stays in main process)
    reader = igzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    # start a process pool
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(polyt_min_length, adapter, barcode_len, umi_len, search_range,
                                                min_qual, whitelist or set(), upstream_region_tsv is not None))

    optional_upstream_cm = open(upstream_region_tsv, 'w') if upstream_region_tsv is not None else nullcontext(None)

    with (open(putative_barcode_file, 'w') as putative_set_file, optional_upstream_cm as upstream_region_file, executor_cm as executor):
        if upstream_region_file is not None:
            upstream_region_file.write('read_id\tsearch_start_to_polyt\n')

        pending = []
        batch = []
        with reader as file:
            for seq_title, seq, qual in FastqGeneralIterator(file):
                seq_title = seq_title.split()[0]
                batch.append((seq_title, seq, qual))
                if len(batch) >= batch_size:
                    pending.append(executor.submit(_worker_process_batch, batch))
                    batch = []

                # keep pending futures bounded and flush completed ones
                if len(pending) >= max_workers * 3:
                    done_future = next(as_completed(pending))
                    pending.remove(done_future)
                    res = _safe_future_result(done_future)
                    result_dict = res

                    for line in result_dict['local_putative_lines']:
                        putative_set_file.write(line)
                    if upstream_region_file is not None:
                        for line in result_dict['local_upstream_region_lines']:
                            upstream_region_file.write(line)
                    for k, v in result_dict['local_putative_freq'].items():
                        putative_barcodes_freq[k] += v
                    for k, v in result_dict['local_anchor_freq'].items():
                        putative_anchors_freq[k] += v
                    polya_found_reads += result_dict['local_polya']
                    reads_with_putative_anchor_count += result_dict['local_anchor']
                    processed += result_dict['local_processed']
                    total_quality_filtered_barcodes += result_dict['local_quality_filtered_barcodes']
                    only_putative_barcodes += result_dict['local_only_putative_reads']
                    both_strands_valid += result_dict['local_both']
                    adapter_found += result_dict['local_adapter']
                    if processed % 500000 == 0:
                        logger.info('Extracted barcodes from %d reads', processed)

            # submit any remaining records
            if batch:
                pending.append(executor.submit(_worker_process_batch, batch))

            # collect remaining futures
            for future in as_completed(pending):
                res = _safe_future_result(future)

                result_dict = res

                for line in result_dict['local_putative_lines']:
                    putative_set_file.write(line)
                if upstream_region_file is not None:
                    for line in result_dict['local_upstream_region_lines']:
                        upstream_region_file.write(line)
                for k, v in result_dict['local_putative_freq'].items():
                    putative_barcodes_freq[k] += v
                for k, v in result_dict['local_anchor_freq'].items():
                    putative_anchors_freq[k] += v
                polya_found_reads += result_dict['local_polya']
                reads_with_putative_anchor_count += result_dict['local_anchor']
                processed += result_dict['local_processed']
                total_quality_filtered_barcodes += result_dict['local_quality_filtered_barcodes']
                only_putative_barcodes += result_dict['local_only_putative_reads']
                both_strands_valid += result_dict['local_both']
                adapter_found += result_dict['local_adapter']
                if processed % 500000 == 0:
                    logger.info('Extracted barcodes from %d reads', processed)

    with open(anchor_freqs_file, 'w') as out:
        for key, value in putative_anchors_freq.items():
            out.write(f'{key},{value}\n')

    prior_count = len(putative_anchors_freq)
    results = _fit_nb_mixture_2component(
        putative_anchors_freq,
        probability_threshold,
        output_csv=probabilities_csv,
        expected_cells=expected_cells,
        max_iter=max_iter,
    )
    cutoff = results['threshold']
    _plot_anchor_probs(putative_anchors_freq, results['probs_real'], results['probs_noise'], cutoff, probability_threshold, probability_figure)
    if cutoff is None:
        cutoff = 0
    putative_anchors_freq = {key: value for key, value in putative_anchors_freq.items() if value >= cutoff}

    # log information
    logger.info('Total reads processed: %d', processed)
    logger.info('polyT+ reads found: %d', polya_found_reads)
    logger.info('Reads with an adapter found: %d', adapter_found)
    logger.info('Total quality filtered barcodes: %d', total_quality_filtered_barcodes)
    logger.info('Passed reads with anchors: %d', reads_with_putative_anchor_count)
    logger.info('Passed reads without anchors: %d', only_putative_barcodes)
    logger.info('Maximum recoverable reads: %d', reads_with_putative_anchor_count + only_putative_barcodes)
    logger.info('Number of putative anchor barcodes before filter: %d', prior_count)
    logger.info('Minimum anchor size calculated from NB is: %d', cutoff)
    logger.info('The number of anchors after filtering is: %d', len(putative_anchors_freq))
    logger.info('Number of unique putative barcodes: %d', len(putative_barcodes_freq))

    putative_anchors_freq = dict(sorted(putative_anchors_freq.items()))
    putative_barcodes_freq = dict(sorted(putative_barcodes_freq.items()))

    return putative_anchors_freq, putative_barcodes_freq
