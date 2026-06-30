from collections import defaultdict
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, as_completed
import lmdb
from isal import igzip 
import subprocess
import os
import logging

logger = logging.getLogger(__name__)

_worker_barcode_to_community = None
_worker_reads_lookup = None
_worker_reads_env = None
_worker_reads_txn = None
_worker_anchor_set = None
_worker_pigz_threads = None


def _lookup_read_sequence(read_id):
    if _worker_reads_txn is not None:
        value = _worker_reads_txn.get(read_id.encode('utf-8'))
        return None if value is None else bytes(value).decode('utf-8')

    return _worker_reads_lookup.get(read_id)

def _worker_process_batch(args):
    chunk, worker_id = args

    local_cell_membership = defaultdict(set)
    local_cell_sizes = defaultdict(int)

    barcode_to_community = _worker_barcode_to_community
    community_to_anchor = _worker_community_to_anchor

    shard = f"assigned_reads.part_{worker_id}.fasta.gz"

    buf = []
    append = buf.append
    flush_threshold = 100000

    with igzip.open(shard, "wt") as out:
        for row in chunk:
            read_id, barcode_string, *_ = row.rstrip("\n").split("\t")

            barcode = None
            community = None

            for bc in barcode_string.split(","):
                community = barcode_to_community.get(bc)
                if community is not None:
                    barcode = bc
                    break

            if barcode is None:
                continue

            local_cell_membership[community].add(read_id)
            local_cell_sizes[community] += 1

            seq = _lookup_read_sequence(read_id)
            if seq is None:
                continue

            append(
                f">{read_id}|barcode={barcode}|community={community}"
                f"|community_anchor={community_to_anchor[community]}\n{seq}"
            )

            if len(buf) >= flush_threshold:
                out.write("\n".join(buf) + "\n")
                buf.clear()

        if buf:
            out.write("\n".join(buf) + "\n")

    return local_cell_membership, local_cell_sizes, shard


def _worker_init(barcode_to_community, community_to_anchor, reads_lookup, anchor_set, pigz_threads):
    global _worker_barcode_to_community
    global _worker_community_to_anchor
    global _worker_reads_lookup
    global _worker_reads_env
    global _worker_reads_txn
    global _worker_anchor_set
    global _worker_pigz_threads

    _worker_barcode_to_community = barcode_to_community
    _worker_community_to_anchor = community_to_anchor
    if isinstance(reads_lookup, str):
        _worker_reads_env = lmdb.open(
            reads_lookup,
            subdir=True,
            readonly=True,
            lock=False,
            readahead=False,
            max_dbs=1,
        )
        _worker_reads_txn = _worker_reads_env.begin()
        _worker_reads_lookup = None
    else:
        _worker_reads_lookup = reads_lookup
        _worker_reads_env = None
        _worker_reads_txn = None
    _worker_anchor_set = anchor_set
    _worker_pigz_threads = pigz_threads


def _submit_batch(executor, chunk, worker_id):
    return executor.submit(_worker_process_batch, (chunk, worker_id))


def assignment(graph_data, reads_lookup, anchor_set, processes=None, batch_size=200000, max_pending=7):
    # initialize variables
    barcode_to_community = defaultdict(str)
    community_to_anchor = defaultdict(list)
    for barcode, community in zip(graph_data['name'], graph_data['community']):
        barcode_to_community[barcode] = community
        if barcode in anchor_set:
            community_to_anchor[community].append(barcode)
    community_to_anchor = {community: ';'.join(anchors) for community, anchors in community_to_anchor.items()}
    cell_membership = defaultdict(set)
    cell_sizes = defaultdict(int)
    used_ids = set()
    shard_files = []
    if processes == None:
        processes = os.cpu_count() - 1 // 2
    else:
        processes = processes // 2
    pending = set()
    worker_id = 0
    logger.info('Starting assignment')
    with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init, initargs=(barcode_to_community, community_to_anchor, reads_lookup, anchor_set, 2)) as executor, open('.putative_barcodes', 'r') as file:
        batch = []
        prior_read = None

        for line in file:
            batch.append(line)
            # make and submit batches
            if len(batch) == batch_size:
                future = _submit_batch(executor, batch, worker_id)
                pending.add(future)
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

        # submit final chunk and process it
        if batch:
            future = _submit_batch(executor, batch, worker_id)
            pending.add(future)
            worker_id += 1

        for fut in as_completed(pending):
            batch_cell_membership, batch_cell_sizes, shard = fut.result()

            shard_files.append(shard)

            for k, v in batch_cell_membership.items():
                cell_membership[k].update(v)
            for k, v in batch_cell_sizes.items():
                cell_sizes[k] += v

    logger.info("Merging FASTA shards")

    subprocess.run(
        "cat assigned_reads.part_*.fasta.gz > assigned_reads.fasta.gz",
        shell=True,
        check=True
    )

    logger.info("Merge complete.")

    # cleanup
    for shard in shard_files:
        os.remove(shard)

    with open('cell_membership.tsv', 'w') as f:
        f.write("community\treads\n")
        for comm, reads in cell_membership.items():
            f.write(f"{comm}\t{' '.join(reads)}\n")

    with open('cell_sizes.tsv', 'w') as f:
        f.write("community\tsize\n")
        for comm, size in cell_sizes.items():
            f.write(f"{comm}\t{size}\n")

    logger.info("Assignment complete.")
