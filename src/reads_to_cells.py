from collections import defaultdict, Counter
from itertools import groupby
import logging
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def assignment(graph_data, putative_frequency):
    """
    Takes the community membership information to assign the reads to cells.

    Args:
        graph_data (pd.DataFrame) : DataFrame of the results of community detection and cell cutoff
        putative_frequency (dict):  Dict of frequency of putative barcodes
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
    multi_assigned = defaultdict(set)

    # track used to see if a read has multiple members
    valid_barcodes_per_read = defaultdict(int)

    with open('cell_sizes.tsv', 'w') as cell_sizes_f, open('cell_membership.tsv', 'w') as cell_membership_f, open(
            '.putative_barcodes', 'r') as putatives_file, open('read_community.tsv', 'w') as read_community_f:
        for read, group in groupby(putatives_file, key=lambda line: line.split('\t', 1)[0]):
            grouped_data = [line.rstrip('\n').split('\t') for line in group]
            valid_barcodes = 0
            communities = []
            barcodes = []
            umis = []
            location_ranks = []
            for row in grouped_data:
                barcode = row[1]
                umi = row[2]
                location_rank = int(row[3])
                if barcode in graph_barcodes:
                    valid_barcodes += 1
                    communities.append(str(barcode_community_lookup[barcode]))
                    barcodes.append(barcode)
                    umis.append(umi)
                    location_ranks.append(location_rank)
            valid_barcodes_per_read[valid_barcodes] += 1
            if valid_barcodes == 1:
                cell_membership[communities[0]].add(read)
                cell_sizes[communities[0]] += 1
                read_community[read].add(str(communities[0]))
            elif valid_barcodes > 1:
                for i in range(valid_barcodes):
                    multi_assigned[read].add((location_ranks[i], communities[i]))

        logger.info('Distribution of valid barcodes per read before rescue: %s', valid_barcodes_per_read)

        # multi-assigned rescue based upon location rank
        all_ranks = []
        used_locations = []
        choices = []
        for read, matches in multi_assigned.items():
            match_dict = {}
            read_ranks = []
            for location_rank, community in matches:
                match_dict[location_rank] = community
                read_ranks.append(location_rank)
                all_ranks.append(location_rank)
            most_likely = sorted(read_ranks, key=lambda x: (abs(x), x))[0]
            choices.append(';'.join([str(x) for x in sorted(read_ranks)]))
            used_locations.append(most_likely)
            community = match_dict[most_likely]
            cell_membership[community].add(read)
            cell_sizes[community] += 1
            read_community[read].add(str(community))
            
        # visualize the locations, maybe try a couple out. Output new shape of valid_barcodes_per_read
        location_counter = Counter(used_locations)

        print(location_counter)
        print(Counter(choices))

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


