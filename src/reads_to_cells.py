from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


def assignment(graph_data):
    """
    Takes the community membership information to assign the reads to cells.

    Args:
        graph_data (pd.DataFrame) : DataFrame of the results of community detection and cell cutoff
    """
    # collect the set of barcodes in the communities for quick search
    graph_barcodes = set(graph_data['name'])
    

    # make defaultdicts for cell membership
    cell_membership = defaultdict(list)
    cell_sizes = defaultdict(int)

    # track used to see if a read has multiple members
    used_reads = set()
    multi_reads = 0

    with open('cell_sizes.csv', 'w') as cell_sizes_f, open('cell_membership.csv', 'w') as cell_membership_f, open('.putative_barcodes', 'r') as putatives_file:
        for line in putatives_file:
            putative_data = line.strip().split('\t')
            read = putative_data[0]
            barcode = putative_data[1]
            if barcode in graph_barcodes:
                cell_sizes[barcode] += 1
                cell_membership[barcode].append(read)
                if read in used_reads:
                    multi_reads += 1
                used_reads.add(read)

        logger.info('Number of multi-assigned reads: %d', multi_reads)
        
        # write out the data
        for key, value in cell_membership.items():
            cell_membership_f.write(f'{key},{" ".join(value)}\n')

        for key, value in cell_sizes.items():
            cell_sizes_f.write(f'{key},{value}\n')
