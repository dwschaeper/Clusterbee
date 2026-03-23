from collections import defaultdict
from itertools import groupby
import logging

logger = logging.getLogger(__name__)


def assignment(graph_data):
    """
    Takes the community membership information to assign the reads to cells.

    Args:
        graph_data (pd.DataFrame) : DataFrame of the results of community detection and cell cutoff
    """
    # relate barcodes to community
    graph_barcodes = graph_data['name'].values.tolist()
    graph_communities = graph_data['community'].tolist()
    barcode_community_lookup = defaultdict(list)
    for barcode, community in zip(graph_barcodes, graph_communities):
        barcode_community_lookup[barcode].append(community)

    # collect the set of barcodes in the communities for quick search
    graph_barcodes = set(graph_barcodes)

    # make defaultdicts for cell membership
    cell_membership = defaultdict(set)
    cell_sizes = defaultdict(int)
    read_community = defaultdict(set)

    # track used to see if a read has multiple members
    valid_barcodes_per_read = defaultdict(int)

    with open('cell_sizes.tsv', 'w') as cell_sizes_f, open('cell_membership.tsv', 'w') as cell_membership_f, open(
            '.putative_barcodes', 'r') as putatives_file, open('read_community.tsv', 'w') as read_community_f:
        for read, group in groupby(putatives_file, key=lambda line: line.split('\t', 1)[0]):
            grouped_data = [line.rstrip('\n').split('\t') for line in group]
            valid_barcodes = 0
            communities = []
            umis = []
            for row in grouped_data:
                barcode = row[1]
                umi = row[2]
                if barcode in graph_barcodes:
                    valid_barcodes += 1
                    communities.append(str(barcode_community_lookup[barcode]))
                    umis.append(umi)
            valid_barcodes_per_read[valid_barcodes] += 1
            if len(communities) == 1:
                cell_membership[communities[0]].add(read)
                cell_sizes[communities[0]] += 1
                read_community[read].add(str(communities[0]))
        
        logger.info('Distribution of valid barcodes per read: %s', valid_barcodes_per_read)

        # write out the data
        cell_membership_f.write('community\treads\n')
        for key, value in cell_membership.items():
            cell_membership_f.write(f'{key}\t{" ".join(value)}\n')

        cell_sizes_f.write('community\tsize\n')
        for key, value in cell_sizes.items():
            cell_sizes_f.write(f'{key}\t{value}\n')

        read_community_f.write('read\tcommunity1\tcommunity2\n')
        for key, value in read_community.items():
            if len(value) > 1:
                comm_string = "\t".join(sorted(list(value)))
                read_community_f.write(f'{key}\t{comm_string}\n')


