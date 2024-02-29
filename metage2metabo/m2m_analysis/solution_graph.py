# Copyright (C) 2019-2024 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade - Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import csv
import json
import logging
import networkx as nx
import os
import sys
import time

from itertools import combinations
from metage2metabo import utils, sbml_management
from metage2metabo.m2m_analysis.taxonomy import get_taxon, extract_taxa
from operator import add

logger = logging.getLogger(__name__)


def graph_analysis(json_file_folder, target_folder_file, output_dir, taxon_file=None, taxonomy_level="phylum"):
    """Run the graph creation on miscoto output

    Args:
        json_file_folder (str): json file or folder containing multiple jsons
        target_folder_file (str): targets file or folder containing multiple sbmls
        output_dir (str): results directory
        taxon_file (str): mpwt taxon file for species in sbml folder
        taxonomy_level (str): taxonomy level, must be: phylum, class, order, family, genus or species.

    Returns:
        str: path to folder containing gml results
    """
    starttime = time.time()
    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#         Solution graph creation             #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    target_paths = utils.file_or_folder(target_folder_file)
    json_paths = utils.file_or_folder(json_file_folder)

    gml_output = os.path.join(output_dir, 'gml')

    if taxon_file:
        taxonomy_output_file = os.path.join(output_dir, 'taxonomy_species.tsv')
        tree_output_file = os.path.join(output_dir, 'taxon_tree.txt')
        extract_taxa(taxon_file, taxonomy_output_file, tree_output_file, taxonomy_level)
    else:
        taxonomy_output_file = None

    create_gml(json_paths, target_paths, output_dir, taxonomy_output_file)

    logger.info(
        "--- Graph runtime %.2f seconds ---\n" % (time.time() - starttime))

    return gml_output


def create_gml(json_paths, target_paths, output_dir, taxon_file=None):
    """Create solution graph from miscoto output and compute stats

    Args:
        json_paths (str): {target: path_to_corresponding_json}
        target_paths (str): {target: path_to_corresponding_sbml}
        output_dir (str): results directory
        taxon_file (str): mpwt taxon file for species in sbml folder
    """
    miscoto_stat_output = os.path.join(output_dir, 'miscoto_stats.txt')
    key_species_stats_output = os.path.join(output_dir,'key_species_stats.tsv')
    key_species_json = os.path.join(output_dir, 'key_species.json')

    gml_output = os.path.join(output_dir, 'gml')

    if not utils.is_valid_dir(gml_output):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    len_min_sol = {}
    len_union = {}
    len_intersection = {}
    len_solution = {}
    len_target = {}

    target_categories = {}
    for target in target_paths:
        target_categories[target] = sbml_management.get_compounds(target_paths[target])

    if taxon_file:
        taxon_named_species, all_taxons = get_taxon(taxon_file)
    else:
        taxon_named_species = None
        all_taxons = None

    key_species_data = {}
    miscoto_stat_output_datas = []

    for target_category in target_categories:
        key_species_data[target_category] = {}
        key_species_data[target_category]['essential_symbionts'] = {}
        key_species_data[target_category]['alternative_symbionts'] = {}
        target_output_gml_path = os.path.join(gml_output, target_category + '.gml')
        with open(json_paths[target_category]) as json_data:
            dicti = json.load(json_data)

        G = nx.Graph()
        added_node = []
        species_weight = {}
        if dicti['still_unprod'] != []:
            logger.error('ERROR ', dicti["still_unprod"], ' is unproducible')
            logger.error('ERROR: Please remove these unproducible targets ({0}) from the targets file, delete output folder ({1}) re-run m2m_analysis'.format(','.join(dicti["still_unprod"]), output_dir))
            sys.exit(1)

        len_target[target_category] = len(dicti['newly_prod']) + len(dicti['still_unprod'])
        len_min_sol[target_category] = len(dicti['bacteria'])
        len_union[target_category] = len(dicti['union_bacteria'])
        len_intersection[target_category] = len(dicti['inter_bacteria'])
        key_species_types = {organism:'ES' if organism in dicti['inter_bacteria'] else 'AS' for organism in dicti['union_bacteria']}
        if taxon_file:
            for taxon in all_taxons:
                key_species_data[target_category]['essential_symbionts'][taxon] = [organism for organism in key_species_types if key_species_types[organism] == 'ES' and taxon_named_species[organism].split('__')[0] == taxon]
                key_species_data[target_category]['alternative_symbionts'][taxon] = [organism for organism in key_species_types if key_species_types[organism] == 'AS' and taxon_named_species[organism].split('__')[0] == taxon]
        else:
            key_species_data[target_category]['essential_symbionts']['data'] = [organism for organism in key_species_types if key_species_types[organism] == 'ES']
            key_species_data[target_category]['alternative_symbionts']['data'] = [organism for organism in key_species_types if key_species_types[organism] == 'AS']

        len_solution[target_category] = len(dicti['enum_bacteria'])
        for sol in dicti['enum_bacteria']:
            if len(dicti['enum_bacteria'][sol]) > 1:
                for species_1, species_2 in combinations(dicti['enum_bacteria'][sol], 2):
                    if species_1 not in added_node:
                        if taxon_file:
                            G.add_node(taxon_named_species[species_1], note=key_species_types[species_1])
                        else:
                            G.add_node(species_1, note=key_species_types[species_1])
                        added_node.append(species_1)
                    if species_2 not in added_node:
                        if taxon_file:
                            G.add_node(taxon_named_species[species_2], note=key_species_types[species_2])
                        else:
                            G.add_node(species_2, note=key_species_types[species_2])
                        added_node.append(species_2)
                    combination_species = '_'.join(sorted([species_1, species_2]))
                    if combination_species not in species_weight:
                        species_weight[combination_species] = 1
                    else:
                        species_weight[combination_species] += 1
                    if taxon_file:
                        G.add_edge(taxon_named_species[species_1], taxon_named_species[species_2], weight=species_weight[combination_species])
                    else:
                        G.add_edge(species_1, species_2, weight=species_weight[combination_species])
            elif len(dicti['enum_bacteria'][sol]) == 1:
                species_1 = dicti['enum_bacteria'][sol][0]
                if species_1 not in added_node:
                    if taxon_file:
                        G.add_node(taxon_named_species[species_1], note=key_species_types[species_1])
                    else:
                        G.add_node(species_1, note=key_species_types[species_1])
                    added_node.append(species_1)

        # Check if all the nodes of G are not isolates.
        if len(G.nodes) == nx.number_of_isolates(G):
            logger.critical(r'/!\ Warning: All the nodes of the solution graph are isolated (they are not connected to other nodes). This lead to powergrasp creating an empty powergraph.')
            logger.critical('So m2m_analysis stops at the solution graph step.')
            sys.exit(1)

        miscoto_stat_output_datas.append([target_category, str(len_target[target_category]), str(len_min_sol[target_category]),
                                str(len_union[target_category]), str(len_intersection[target_category]),
                                str(len_solution[target_category])])
        logger.info('######### Graph of ' + target_category + ' #########')
        logger.info('Number of nodes: ' + str(G.number_of_nodes()))
        logger.info('Number of edges: ' + str(G.number_of_edges()))
        nx.write_gml(G, target_output_gml_path)

    with open(miscoto_stat_output, 'w') as stats_output:
        statswriter = csv.writer(stats_output, delimiter="\t")
        statswriter.writerow(['categories', 'nb_target', 'size_min_sol', 'size_union', 'size_intersection', 'size_enum'])
        for miscoto_stat_output_data in miscoto_stat_output_datas:
            statswriter.writerow(miscoto_stat_output_data)

    with open(key_species_json, 'w') as json_output:
        json.dump(key_species_data, json_output, indent=4)

    with open(key_species_stats_output, 'w') as key_stat_output:
        key_stats_writer = csv.writer(key_stat_output, delimiter='\t')
        if all_taxons:
            key_stats_writer.writerow(['target_categories', 'key_group', *sorted(all_taxons), 'Sum'])
        else:
            key_stats_writer.writerow(['target_categories', 'key_group', 'data', 'Sum'])
        for target in key_species_data:
            if all_taxons:
                essential_counts = [len(key_species_data[target]['essential_symbionts'][taxon]) for taxon in sorted(all_taxons)]
                alternative_counts = [len(key_species_data[target]['alternative_symbionts'][taxon]) for taxon in sorted(all_taxons)]
            else:
                essential_counts = [len(key_species_data[target]['essential_symbionts']['data'])]
                alternative_counts = [len(key_species_data[target]['alternative_symbionts']['data'])]
            key_counts = list(map(add, essential_counts, alternative_counts))
            key_stats_writer.writerow([target, 'key_species', *key_counts, sum(key_counts)])
            key_stats_writer.writerow([target, 'essential_symbionts', *essential_counts, sum(essential_counts)])
            key_stats_writer.writerow([target, 'alternative_symbionts', *alternative_counts, sum(alternative_counts)])