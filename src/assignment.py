from collections import defaultdict
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import os
import logging

logger = logging.getLogger(__name__)

_worker_barcode_to_community = None
_worker_valid_barcodes = None
_worker_reads_lookup = None
_worker_pigz_threads = None


def _write_fasta_shard(read_assignment, worker_id):
    outfile = f"assigned_reads.part_{worker_id}.fasta.gz"

    proc = subprocess.Popen(
        ["pigz", "-p", str(_worker_pigz_threads), "-c"],
        stdin=subprocess.PIPE,
        stdout=open(outfile, "wb"),
        text=True
    )

    buf = []
    append = buf.append
    flush_threshold = 100000

    for read_id, (barcode, community) in read_assignment.items():
        seq = _worker_reads_lookup.get(read_id)
        if seq is None:
            continue

        append(f">{read_id}|barcode={barcode}|community={community}\n{seq}")

        if len(buf) >= flush_threshold:
            proc.stdin.write("\n".join(buf) + "\n")
            buf.clear()

    if buf:
        proc.stdin.write("\n".join(buf) + "\n")

    proc.stdin.close()
    proc.wait()

    return outfile


def _worker_process_batch(args):
    chunk, worker_id = args

    local_cell_membership = defaultdict(set)
    local_cell_sizes = defaultdict(int)
    read_assignment = {}

    for read_id, group in groupby(chunk, key=lambda x: x.split('\t', 1)[0]):
        rows = [line.rstrip('\n').split('\t') for line in group]

        matches = []
        for row in rows:
            barcode = row[1]
            location_rank = int(row[3])

            if barcode in _worker_valid_barcodes:
                community = _worker_barcode_to_community[barcode]
                matches.append((location_rank, barcode, community))

        # if only one match, assign
        if len(matches) == 1:
            _, barcode, community = matches[0]

        # if multiple, take best position
        elif len(matches) > 1:
            location_rank, barcode, community = sorted(matches, key=lambda x: (abs(x[0]), x[0]))[0]
        else:
            continue

        local_cell_membership[community].add(read_id)
        local_cell_sizes[community] += 1
        read_assignment[read_id] = (barcode, community)

    shard = _write_fasta_shard(read_assignment, worker_id)

    return local_cell_membership, local_cell_sizes, shard


def _worker_init(barcode_to_community, valid_barcodes, reads_lookup, pigz_threads):
    global _worker_barcode_to_community
    global _worker_valid_barcodes
    global _worker_reads_lookup
    global _worker_pigz_threads

    _worker_barcode_to_community = barcode_to_community
    _worker_valid_barcodes = valid_barcodes
    _worker_reads_lookup = reads_lookup
    _worker_pigz_threads = pigz_threads


def _submit_batch(executor, chunk, worker_id):
    return executor.submit(_worker_process_batch, (chunk, worker_id))


def assignment(graph_data, reads_lookup, processes=None, batch_size=10000, max_pending=7):
    # initialize variables
    barcode_to_community = dict(zip(graph_data['name'], graph_data['community']))
    valid_barcodes = set(barcode_to_community.keys())
    cell_membership = defaultdict(set)
    cell_sizes = defaultdict(int)
    shard_files = []
    if processes == None:
        processes = os.cpu_count() - 1
    pending = set()
    worker_id = 0
    logger.info('Starting assignment')
    with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init, initargs=(barcode_to_community, valid_barcodes, reads_lookup, 1)) as executor, open('.putative_barcodes', 'r') as file:
        batch = []
        current_read = None

        for line in file:
            read_id = line.split('\t', 1)[0]

            if current_read is None:
                current_read = read_id

            # make and submit batches
            if read_id != current_read and len(batch) >= batch_size:  # ensure the putatives for all reads are within a single batch
                future = _submit_batch(executor, batch, worker_id)
                pending.add(future)
                worker_id += 1

                batch = []
                current_read = read_id

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

            batch.append(line)

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
