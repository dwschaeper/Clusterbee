from collections import defaultdict
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

    # check to see that a barcode only has one associated community
    multi = 0
    for key, value in barcode_community_lookup.items():
        if len(value) > 1:
            multi += 1
    logger.info('Number of barcodes that belong to multiple communities: %d', multi)

    # collect the set of barcodes in the communities for quick search
    graph_barcodes = set(graph_barcodes)
    

    # make defaultdicts for cell membership
    cell_membership = defaultdict(set)
    cell_sizes = defaultdict(int)
    read_community = defaultdict(set)

    # track used to see if a read has multiple members
    used_reads = set()
    multi_reads = 0

    with open('cell_sizes.tsv', 'w') as cell_sizes_f, open('cell_membership.tsv', 'w') as cell_membership_f, open('.putative_barcodes', 'r') as putatives_file, open('read_community.tsv', 'w') as read_community_f:
        for line in putatives_file:
            putative_data = line.strip().split('\t')
            read = putative_data[0]
            barcode = putative_data[1]
            if barcode in graph_barcodes:
                community = barcode_community_lookup[barcode]
                # track only once per community
                if read not in cell_membership[community[0]]:
                    cell_membership[community[0]].add(read)
                    cell_sizes[community[0]] += 1
                    read_community[read].add(str(community[0]))
                    # check to see if read is present in another community
                    if read in used_reads:
                        multi_reads += 1
                    used_reads.add(read)

        logger.info('Number of multi-assigned reads: %d', multi_reads)
        
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


