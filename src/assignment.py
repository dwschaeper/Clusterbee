from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from isal import igzip 
import subprocess
import os
import logging

logger = logging.getLogger(__name__)

_worker_barcode_to_community = None
_worker_community_to_anchor = None
_worker_anchor_set = None

def _worker_process_batch(args):
    batch, worker_id = args

    local_cell_membership = defaultdict(set)
    local_cell_sizes = defaultdict(int)
    local_position_dist = defaultdict(int)

    barcode_to_community = _worker_barcode_to_community
    community_to_anchor = _worker_community_to_anchor

    shard = f"assigned_reads.part_{worker_id}.fastq.gz"

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
            local_position_dist[pos] += 1

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

    return local_cell_membership, local_cell_sizes, local_position_dist, shard


def _worker_init(barcode_to_community, community_to_anchor, anchor_set):
    global _worker_barcode_to_community
    global _worker_community_to_anchor
    global _worker_anchor_set

    _worker_barcode_to_community = barcode_to_community
    _worker_community_to_anchor = community_to_anchor
    _worker_anchor_set = anchor_set


def _submit_batch(executor, batch, worker_id):
    return executor.submit(_worker_process_batch, (batch, worker_id))


def assignment(graph_data, anchor_set, processes=None, batch_size=250000, max_pending=15):
    # initialize variables
    barcode_to_community = defaultdict(str)
    community_to_anchor = defaultdict(list)
    position_dict = defaultdict(int)
    for barcode, community in zip(graph_data['name'], graph_data['community']):
        barcode_to_community[barcode] = community
        if barcode in anchor_set:
            community_to_anchor[community].append(barcode)
    community_to_anchor = {community: ';'.join(anchors) for community, anchors in community_to_anchor.items()}
    cell_membership = defaultdict(set)
    cell_sizes = defaultdict(int)
    shard_files = []
    if processes == None:
        processes = (os.cpu_count() - 1)
    else:
        processes = processes
    pending = set()
    worker_id = 0
    logger.info('Starting assignment')
    with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init, initargs=(barcode_to_community, community_to_anchor, anchor_set)) as executor, open('.putative_barcodes', 'r') as file:
        batch = []

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
                        batch_cell_membership, batch_cell_sizes, batch_position_dict, shard = fut.result()

                        shard_files.append(shard)

                        for k, v in batch_cell_membership.items():
                            cell_membership[k].update(v)
                        for k, v in batch_cell_sizes.items():
                            cell_sizes[k] += v
                        for k, v in batch_position_dict.items():
                            position_dict[k] += v

                        completed.append(fut)

                        if len(pending) - len(completed) < max_pending:
                            break

                    for fut in completed:
                        pending.remove(fut)

        # submit final batch and process it
        if batch:
            future = _submit_batch(executor, batch, worker_id)
            pending.add(future)
            worker_id += 1

        for fut in as_completed(pending):
            batch_cell_membership, batch_cell_sizes, batch_position_dict, shard = fut.result()

            shard_files.append(shard)

            for k, v in batch_cell_membership.items():
                cell_membership[k].update(v)
            for k, v in batch_cell_sizes.items():
                cell_sizes[k] += v

    logger.info("Merging FASTA shards")

    subprocess.run(
        "cat assigned_reads.part_*.fastq.gz > assigned_reads.fastq.gz",
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
    logger.info(position_dict)
