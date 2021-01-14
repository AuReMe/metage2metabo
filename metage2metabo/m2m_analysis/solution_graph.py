import csv
import json
import logging
import networkx as nx
import os
import sys
import time

from itertools import combinations
from metage2metabo import utils, sbml_management
from metage2metabo.m2m_analysis.taxonomy import detect_taxon_species, get_taxon, extract_taxa

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


def detect_key_species(json_elements, all_taxons, taxon_named_species=None):
    """Detect key species (essential and alternative symbionts) from the miscoto json

    Args:
        json_elements (dict): miscoto results in a json dictionary
        all_taxons (list): all taxon in the dataset
        taxon_named_species (dict): {species_ID: species_named_after_taxon}

    Returns:
        key_species (dict): {taxon: [species_1, species_2]}
        essential_symbionts (dict): {taxon: [species_1, species_2]}
        alternative_symbionts (dict): {taxon: [species_1, species_2]}
    """
    if taxon_named_species:
        unions = {taxon_named_species[species_union]: species_union
            		for species_union in json_elements["union_bacteria"]}
        intersections = {taxon_named_species[species_intersection]: species_intersection
            		for species_intersection in json_elements["inter_bacteria"]}
    else:
        unions = json_elements["union_bacteria"]
        intersections = json_elements["inter_bacteria"]

    if taxon_named_species:
        key_species = detect_taxon_species(unions, all_taxons)
        essential_symbionts = detect_taxon_species(intersections, all_taxons)

        alternative_symbionts = {}
        for taxon in key_species:
            if taxon in essential_symbionts:
                alternative_symbionts[taxon] = list(set(key_species[taxon]) - set(essential_symbionts[taxon]))
            else:
                alternative_symbionts[taxon] = list(set(key_species[taxon]))
        for taxon in all_taxons:
            if taxon not in alternative_symbionts:
                alternative_symbionts[taxon] = []
    else:
        key_species = {}
        essential_symbionts = {}
        alternative_symbionts = {}
        key_species["data"] = unions
        essential_symbionts["data"] = intersections
        alternative_symbionts["data"] = list(set(unions) - set(intersections))

    return key_species, essential_symbionts, alternative_symbionts


def create_stat_species(target_category, json_elements, key_stats_writer, key_sup_writer, taxon_named_species=None, all_taxons=None):
    """Write stats on key species (essential and alternative symbionts) from the miscoto json

    Args:
        target_category (str): name of a target file (without extension)
        json_elements (dict): miscoto results in a json dictionary
        key_stats_writer (csv.writer): writer for stats of each group (key species) in each taxon
        key_sup_writer (csv.writer): writer for all key species in each taxon
        taxon_named_species (dict): {species_ID: species_named_after_taxon}
        all_taxons (list): all taxon in the dataset
    """
    key_stone_species, essential_symbionts, alternative_symbionts = detect_key_species(json_elements, all_taxons, taxon_named_species)

    key_stone_counts = [len(key_stone_species[organism]) for organism in sorted(list(key_stone_species.keys()))]
    key_stats_writer.writerow([target_category, "key_species"] + key_stone_counts + [sum(key_stone_counts)])

    essential_symbiont_counts = [len(essential_symbionts[organism]) for organism in sorted(list(essential_symbionts.keys()))]
    key_stats_writer.writerow([target_category, "essential_symbionts"] + essential_symbiont_counts + [sum(essential_symbiont_counts)])

    alternative_symbiont_counts = [len(alternative_symbionts[organism]) for organism in sorted(list(alternative_symbionts.keys()))]
    key_stats_writer.writerow([target_category, "alternative_symbionts"] + alternative_symbiont_counts + [sum(alternative_symbiont_counts)])

    if all_taxons:
        for taxon in sorted(all_taxons):
            key_sup_writer.writerow([target_category, "key_stone_species", taxon] + key_stone_species[taxon])
            if taxon in essential_symbionts:
                key_sup_writer.writerow([target_category, "essential_symbionts", taxon] + essential_symbionts[taxon])
            else:
                key_sup_writer.writerow([target_category, "essential_symbionts", taxon])
            if taxon in alternative_symbionts:
                key_sup_writer.writerow([target_category, "alternative_symbionts", taxon] + alternative_symbionts[taxon])
            else:
                key_sup_writer.writerow([target_category, "alternative_symbionts", taxon])
    else:
        key_sup_writer.writerow([target_category, "key_stone_species", "data"] + key_stone_species["data"])
        key_sup_writer.writerow([target_category, "essential_symbionts", "data"] + essential_symbionts["data"])
        key_sup_writer.writerow([target_category, "alternative_symbionts", "data"] + alternative_symbionts["data"])


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
    key_species_supdata_output = os.path.join(output_dir, 'key_species_supdata.tsv')

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

    miscoto_stat_output_datas = []
    with open(key_species_stats_output, 'w') as key_stats_file, open(key_species_supdata_output, 'w') as key_sup_file:
        key_stats_writer = csv.writer(key_stats_file, delimiter='\t')
        if all_taxons:
            key_stats_writer.writerow(['target_categories', 'key_group', *sorted(all_taxons), 'Sum'])
        else:
            key_stats_writer.writerow(['target_categories', 'key_group', 'data', 'Sum'])
        key_sup_writer = csv.writer(key_sup_file, delimiter='\t')
        for target_category in target_categories:
            target_output_gml_path = os.path.join(gml_output, target_category + '.gml')
            with open(json_paths[target_category]) as json_data:
                dicti = json.load(json_data)
            create_stat_species(target_category, dicti, key_stats_writer, key_sup_writer, taxon_named_species, all_taxons)
            G = nx.Graph()
            added_node = []
            species_weight = {}
            if dicti['still_unprod'] != []:
                logger.warning('ERROR ', dicti["still_unprod"], ' is unproducible')
            len_target[target_category] = len(dicti['newly_prod']) + len(dicti['still_unprod'])
            len_min_sol[target_category] = len(dicti['bacteria'])
            len_union[target_category] = len(dicti['union_bacteria'])
            len_intersection[target_category] = len(dicti['inter_bacteria'])
            key_species_types = {organism:'ES' if organism in dicti['inter_bacteria'] else 'AS' for organism in dicti['union_bacteria']}
            len_solution[target_category] = len(dicti['enum_bacteria'])
            for sol in dicti['enum_bacteria']:
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
