"""
Assign putative reads to graph communities and emit per-cell outputs.

This module parallelizes read-to-community assignment using the putative barcodes
written during barcode extraction.
"""

from collections import defaultdict
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
import logging
import os
import subprocess
from typing import Any

from isal import igzip

import config

logger = logging.getLogger(__name__)

_worker_barcode_to_community = None
_worker_community_to_anchor = None
_worker_anchor_set = None


def _worker_process_batch(args: tuple[list[str], int]) -> tuple[dict[int, set[str]], dict[int, int], dict[int, int], str]:
    """
    Process a batch of putative-read records in a worker process.

    Args:
        args (tuple[list[str], int]): Tuple of ``(batch_lines, worker_id)``
            where each line is
            ``read_id\tbarcode_candidates\tumi_candidates\tread_seq\tread_qual``.

    Returns:
        tuple[dict[int, set[str]], dict[int, int], dict[int, int], str]:
            Aggregated cell membership, cell sizes, barcode-position
            distribution, and output shard path.
    """
    batch, worker_id = args

    local_cell_membership = defaultdict(set)
    local_cell_sizes = defaultdict(int)

    barcode_to_community = _worker_barcode_to_community
    community_to_anchor = _worker_community_to_anchor

    shard = os.path.join(config.GLOBAL_OUTDIR, f".assigned_reads.part_{worker_id}.fastq.gz")

    buf = []
    append = buf.append
    flush_threshold = 100000

    with igzip.open(shard, "wt") as out:
        for row in batch:
            read_id, barcode_string, umi_string, seq, seq_qual  = row.rstrip("\n").split("\t")

            barcode = None
            community = None
            
            pos = 0
            for bc in barcode_string.split(","):
                community = barcode_to_community.get(bc)
                if community is not None:
                    barcode = bc
                    break
                pos += 1

            if barcode is None:
                continue

            local_cell_membership[community].add(read_id)
            local_cell_sizes[community] += 1

            append(
                f"@{read_id}|barcode={barcode}|community={community}|community_anchor={community_to_anchor[community]}\n"
                f"{seq}\n"
                f"+\n"
                f"{seq_qual}"
            )

            if len(buf) >= flush_threshold:
                out.write("\n".join(buf) + "\n")
                buf.clear()

        if buf:
            out.write("\n".join(buf) + "\n")

    return local_cell_membership, local_cell_sizes, shard


def _worker_init(barcode_to_community: dict[str, int], community_to_anchor: dict[int, str], anchor_set: set[str]) -> None:
    """
    Initialize immutable lookup tables for each worker process.

    Args:
        barcode_to_community (dict[str, int]): Mapping from barcode sequence
            to community id.
        community_to_anchor (dict[int, str]): Mapping from community id to
            joined anchor list.
        anchor_set (set[str]): Set of anchors retained after filtering.
    """
    global _worker_barcode_to_community
    global _worker_community_to_anchor
    global _worker_anchor_set

    _worker_barcode_to_community = barcode_to_community
    _worker_community_to_anchor = community_to_anchor
    _worker_anchor_set = anchor_set


def _submit_batch(executor: ProcessPoolExecutor, batch: list[str], worker_id: int) -> Future[tuple[dict[int, set[str]], dict[int, int], dict[int, int], str]]:
    """
    Submit a worker batch to the process pool.

    Args:
        executor (ProcessPoolExecutor): Active process pool.
        batch (list[str]): Input batch lines.
        worker_id (int): Worker shard id for output naming.

    Returns:
        Future[tuple[dict[int, set[str]], dict[int, int], dict[int, int], str]]:
            Future yielding batch aggregates and shard path.
    """
    return executor.submit(_worker_process_batch, (batch, worker_id))


def assignment(graph_data: dict[str, list[Any]], anchor_set: set[str], processes: int | None = None,
              batch_size: int = 250000, max_pending: int = 15) -> tuple[dict[int, set[str]], dict[int, int], dict[int, int]]:
    """
    Assign putative reads to graph communities and write assignment outputs.

    Args:
        graph_data (dict[str, list[Any]]): Graph table with at least ``name``
            and ``community`` arrays.
        anchor_set (set[str]): Whitelist of surviving anchor barcodes.
        processes (int | None): Number of worker processes. Defaults to
            ``cpu_count() - 1``.
        batch_size (int): Number of putative-read lines per submitted task.
        max_pending (int): Maximum in-flight futures before partial drain.

    Returns:
        tuple[dict[int, set[str]], dict[int, int], dict[int, int]]:
            ``(cell_membership, cell_sizes, barcode_position_distribution)``.
    """
    logger.section("ASSIGNING READS TO CELLS")
    logger.info("Starting assignment")

    barcode_to_community: dict[str, int] = {}
    community_to_anchor: dict[int, list[str]] = defaultdict(list)
    for barcode, community in zip(graph_data["name"], graph_data["community"]):
        barcode_to_community[barcode] = int(community)
        if barcode in anchor_set:
            community_to_anchor[int(community)].append(barcode)

    community_to_anchor_joined: dict[int, str] = {
        community: ";".join(anchors) for community, anchors in community_to_anchor.items()
    }
    cell_membership: dict[int, set[str]] = defaultdict(set)
    cell_sizes: dict[int, int] = defaultdict(int)
    shard_files: list[str] = []
    if processes is None:
        processes = max(1, (os.cpu_count() or 2) - 1)

    pending: set[Future[tuple[dict[int, set[str]], dict[int, int], dict[int, int], str]]] = set()
    worker_id = 0

    putative_barcode_file = os.path.join(config.GLOBAL_OUTDIR, ".putative_barcodes")
    fastq_file = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_assigned_reads.fastq.gz')
    cell_membership_tsv = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_cell_membership.tsv')
    cell_sizes_tsv = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_cell_sizes.tsv')

    with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init, 
        initargs=(barcode_to_community, community_to_anchor_joined, anchor_set)) as executor, open(putative_barcode_file, "r") as file:
        batch = []

        for line in file:
            batch.append(line)
            if len(batch) == batch_size:
                pending.add(_submit_batch(executor, batch, worker_id))
                worker_id += 1
                batch = []

                if len(pending) >= max_pending:
                    completed = []
                    for fut in as_completed(pending):
                        batch_cell_membership, batch_cell_sizes, shard = fut.result()
                        shard_files.append(shard)

                        for k, v in batch_cell_membership.items():
                            cell_membership[k].update(v)
                        for k, v in batch_cell_sizes.items():
                            cell_sizes[k] += v

                        completed.append(fut)

                        if len(pending) - len(completed) < max_pending:
                            break

                    for fut in completed:
                        pending.remove(fut)

        # submit final batch and process it
        if batch:
            pending.add(_submit_batch(executor, batch, worker_id))

        for fut in as_completed(pending):
            batch_cell_membership, batch_cell_sizes, shard = fut.result()
            shard_files.append(shard)

            for k, v in batch_cell_membership.items():
                cell_membership[k].update(v)
            for k, v in batch_cell_sizes.items():
                cell_sizes[k] += v

    logger.info("Merging FASTA shards")

    if config.GLOBAL_OUTDIR[-1] == '/':
        command = f"cat {config.GLOBAL_OUTDIR}.assigned_reads.part_*.fastq.gz > {fastq_file}"
    else:
        command = f"cat {config.GLOBAL_OUTDIR}/.assigned_reads.part_*.fastq.gz > {fastq_file}"

    subprocess.run(
        command,
        shell=True,
        check=True
    )

    logger.info("Merge complete.")

    # cleanup
    for shard in shard_files:
        os.remove(shard)

    with open(cell_membership_tsv, 'w') as f:
        f.write("community\treads\n")
        for comm, reads in cell_membership.items():
            f.write(f"{comm}\t{' '.join(reads)}\n")

    with open(cell_sizes_tsv, 'w') as f:
        f.write("community\tsize\n")
        for comm, size in cell_sizes.items():
            f.write(f"{comm}\t{size}\n")

    logger.info("Assignment complete.")
