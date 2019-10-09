#!/usr/bin/env python
# Copyright (c) 2019, Clemence Frioux <clemence.frioux@inria.fr>
#
# This file is part of metage2metabo.
#
# metage2metabo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# metage2metabo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with metage2metabo.  If not, see <http://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import miscoto
import networkx as nx
import os
import padmet
import powergrasp
import sys
import time
import subprocess
import logging

from ete3 import NCBITaxa
from itertools import combinations
from metage2metabo import utils, sbml_management

recipes = '''{
	"clingo multithreading": 4
}'''
powergrasp.recipe.Recipe.from_lines(recipes)

logger = logging.getLogger(__name__)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)

def file_or_folder(variable_folder_file):
	file_folder_paths = {}

	if os.path.isfile(variable_folder_file):
		filename = os.path.splitext(os.path.basename(variable_folder_file))[0]
		file_folder_paths[filename] = variable_folder_file

	elif os.path.isdir(variable_folder_file):
		for file_from_folder in os.listdir(variable_folder_file):
			filename = os.path.splitext(os.path.basename(file_from_folder))[0]
			file_folder_paths[filename] = variable_folder_file + '/' + file_from_folder

	return file_folder_paths

def run_analysis_workflow(sbml_folder, target_folder_file, seed_file, output_dir, taxon_file, oog_jar, host_file=None):
	enumeration_analysis(seed_file, sbml_folder, target_folder_file, output_dir, host_file)

	for miscoto_json in miscoto_jsons:
		miscoto_analysis = miscoto_json[1]
		miscoto_stat_output = miscoto_analysis + '/' + 'miscoto_stats.txt'
		gml_output = miscoto_analysis + '/' + 'graph.gml'
		bbl_output = miscoto_analysis + '/' + 'graph.bbl'
		svg_output = miscoto_analysis + '/' + 'graph.svg'

		create_gml(target_file, miscoto_json[0], miscoto_analysis, taxon_file)

	powergraph_analysis(output_dir + '/gml', output_dir, oog_jar)


def enumeration(seed_file, sbml_folder, target_file, output_json, host_file):
	results = miscoto.run_mincom(option='soup', bacteria_dir=sbml_folder,
						targets_file=target_file, seeds_file=seed_file, host_file=host_file,
                		intersection=True, enumeration=True, union=True, optsol=True)

	miscoto.utils.results_to_json(results, output_json)

	return output_json


def enumeration_analysis(seed_file, sbml_folder, target_file, output_dir, host_file):
	target_paths = file_or_folder(target_file)

	miscoto_jsons = {}
	for target_path in target_paths:
		target_pathname = target_path[target_path]
		output_json = output_dir + '/' + 'json' + '/' + target_path + '.json'
		if not utils.is_valid_dir(output_json):
			logger.critical("Impossible to access/create output directory")
			sys.exit(1)
		miscoto_json = enumeration(seed_file, sbml_folder, target_pathname, output_json, host_file)
		miscoto_jsons[target_path] = miscoto_json

	return miscoto_jsons


def stat_analysis(json_file_folder, output_dir, taxon_file=None):
	miscoto_stat_output = output_dir + '/' + 'miscoto_stats.txt'
	key_species_stats_output = output_dir + '/' + 'key_species_stats.tsv'
	key_species_supdata_output = output_dir + '/' + 'key_species_supdata.tsv'
	json_paths = file_or_folder(json_file_folder)

	if taxon_file:
		phylum_output_file = output_dir + '/taxon_phylum.tsv'
		tree_output_file = output_dir + '/taxon_tree.txt'
		extract_taxa(taxon_file, phylum_output_file, tree_output_file)
		phylum_species, all_phylums = get_phylum(phylum_output_file)

	with open(key_species_stats_output, "w") as key_stats_file:
		key_stats_writer = csv.writer(key_stats_file, delimiter='\t')
		if all_phylums:
			key_stats_writer.writerow(['target_categories', 'key_stones_group', *sorted(all_phylums), 'Sum'])
		else:
			key_stats_writer.writerow(['target_categories', 'key_stones_group', 'data', 'Sum'])
		with open(key_species_supdata_output,"w") as key_sup_file:
			key_sup_writer = csv.writer(key_sup_file, delimiter='\t')

			for json_path in json_paths:
				with open(json_paths[json_path]) as json_data:
					json_elements = json.load(json_data)
				create_stat_species(json_path, json_elements, key_stats_writer, key_sup_writer, phylum_species, all_phylums)

				with open(miscoto_stat_output,"w") as stats_output:
					statswriter = csv.writer(stats_output, delimiter='\t')
					statswriter.writerow(['categories', 'nb_target', 'size_min_sol', 'size_union', 'size_intersection', 'size_enum'])
					statswriter.writerow([json_path, str(len(json_elements["newly_prod"]) + len(json_elements["still_unprod"])), str(len(json_elements["bacteria"])),
										str(len(json_elements["union_bacteria"])), str(len(json_elements["inter_bacteria"])), str(len(json_elements["enum_bacteria"]))])

def graph_analysis(json_file_folder, output_folder, target_file, taxon_file):
	json_paths = file_or_folder(json_file_folder)

	create_gml(target_file, json_paths, output_folder, taxon_file)

def powergraph_analysis(graph_input, output_folder, oog_jar):
	gml_paths = file_or_folder(graph_input)

	for gml_path in gml_paths:
		bbl_output = output_folder + '/' + 'bbl' + '/' + gml_path + '.bbl'
		svg_output = output_folder + '/' + 'svg' + '/' + gml_path + '.svg'
		graph_input = gml_paths[gml_path]
		compression(graph_input, bbl_output)
		bbl_to_svg(oog_jar, bbl_output, svg_output)

def extract_taxa(mpwt_taxon_file, taxon_output_file, tree_output_file):
	ncbi = NCBITaxa()

	taxon_ids = []

	phylum_count = {}
	with open(taxon_output_file, 'w') as phylum_file:
		csvwriter = csv.writer(phylum_file, delimiter='\t')
		csvwriter.writerow(['species', 'taxid', 'phylum_number', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
		with open(mpwt_taxon_file, 'r') as taxon_file:
			csvfile = csv.reader(taxon_file, delimiter='\t')
			for line in csvfile:
				if 'taxon' not in line[1]:
					taxon_ids.append(line[1])
					lineage = ncbi.get_lineage(line[1])
					lineage2ranks = ncbi.get_rank(lineage)
					names = ncbi.get_taxid_translator(lineage)
					ranks2lineage = dict((rank, names[taxid]) for (taxid, rank) in lineage2ranks.items())
					ranks = [ranks2lineage.get(rank, 'no_information') for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']]
					if ranks[1] != 'no_information':
						phylum = ranks[1][:4]
					else:
						phylum = 'no_information'
					if phylum not in phylum_count :
						phylum_count[phylum] = 1
					elif phylum == 'no_information':
						phylum_count[phylum] = ''
					else:
						phylum_count[phylum] += 1
					row = [line[0], line[1]] + [phylum + str(phylum_count[phylum])] + ranks
					csvwriter.writerow(row)

	tree = ncbi.get_topology(taxon_ids)

	with open(tree_output_file, 'w') as tree_file:
		tree_file.write(tree.get_ascii(attributes=["sci_name", "rank"]))


def detect_phylum_key_species(phylums, all_phylums):
	count_phylumns = {}
	for phylum in phylums:
		if phylum[:4] not in count_phylumns:
			count_phylumns[phylum[:4]] = [phylums[phylum]]
		else:
			count_phylumns[phylum[:4]].append(phylums[phylum])

	for phylum in all_phylums:
		if phylum not in count_phylumns:
			count_phylumns[phylum[:4]] = []

	return count_phylumns


def detect_key_species(json_elements, all_phylums, phylum_species=None):
	if phylum_species:
		unions = {phylum_species[species_union]: species_union for species_union in json_elements["union_bacteria"]}
		intersections = {phylum_species[species_intersection]: species_intersection for species_intersection in json_elements["inter_bacteria"]}
	else:
		unions = json_elements["union_bacteria"]
		intersections = json_elements["inter_bacteria"]

	if phylum_species:
		key_stone_species = detect_phylum_key_species(unions, all_phylums)
		essential_symbionts = detect_phylum_key_species(intersections, all_phylums)
		alternative_symbionts = {}
		for phylum in key_stone_species:
			alternative_symbionts[phylum] = list(set(key_stone_species[phylum]) - set(essential_symbionts[phylum]))
		for phylum in all_phylums:
			if phylum not in alternative_symbionts:
				alternative_symbionts[phylum] = []
	else:
		key_stone_species = {}
		essential_symbionts = {}
		alternative_symbionts = {}
		key_stone_species['data'] = unions
		essential_symbionts['data'] = intersections
		alternative_symbionts['data'] = list(set(unions) - set(intersections))

	return key_stone_species, essential_symbionts, alternative_symbionts


def create_stat_species(target_categories, json_elements, key_stats_writer, key_sup_writer, phylum_species=None, all_phylums=None):
	key_stone_species, essential_symbionts, alternative_symbionts = detect_key_species(json_elements, all_phylums, phylum_species)

	key_stone_counts = [len(key_stone_species[phylum]) for phylum in sorted(list(key_stone_species.keys()))]
	key_stats_writer.writerow([target_categories, 'key_stone_species'] + key_stone_counts + [sum(key_stone_counts)])

	essential_symbiont_counts = [len(essential_symbionts[phylum]) for phylum in sorted(list(essential_symbionts.keys()))]
	key_stats_writer.writerow([target_categories, 'essential_symbionts'] + essential_symbiont_counts + [sum(essential_symbiont_counts)])

	alternative_symbiont_counts = [len(alternative_symbionts[phylum]) for phylum in sorted(list(alternative_symbionts.keys()))]
	key_stats_writer.writerow([target_categories, 'alternative_symbionts'] + alternative_symbiont_counts + [sum(alternative_symbiont_counts)])

	if all_phylums:
		for phylum in sorted(all_phylums):
				key_sup_writer.writerow([target_categories, 'key_stone_species', phylum] + key_stone_species[phylum])
				key_sup_writer.writerow([target_categories, 'essential_symbionts', phylum] + essential_symbionts[phylum])
				key_sup_writer.writerow([target_categories, 'alternative_symbionts', phylum] + alternative_symbionts[phylum])
	else:
		key_sup_writer.writerow([target_categories, 'key_stone_species', 'data'] + key_stone_species[phylum])
		key_sup_writer.writerow([target_categories, 'essential_symbionts', 'data'] + essential_symbionts[phylum])
		key_sup_writer.writerow([target_categories, 'alternative_symbionts', 'data'] + alternative_symbionts[phylum])

def get_phylum(phylum_file):
	phylum_species = {}
	all_phylums = []
	with open(phylum_file, 'r') as phylum_file:
		phylum_reader = csv.reader(phylum_file, delimiter='\t', quotechar='|')
		for row in phylum_reader:
			phylum_species[row[0]] = row[2]
			if row[2][:4] not in all_phylums:
				if 'no_information' not in row[2] and 'phylum_number' not in row[2]:
					all_phylums.append(row[2][:4])

	return phylum_species, all_phylums

def create_gml(targets, jsons, output_dir, taxon_file=None):
	miscoto_stat_output = output_dir + '/' + 'miscoto_stats.txt'
	key_species_stats_output = output_dir + '/' + 'key_species_stats.tsv'
	key_species_supdata_output = output_dir + '/' + 'key_species_supdata.tsv'

	gml_output = output_dir + '/' + 'gml' + '/'

	len_min_sol = {}
	len_union = {}
	len_intersection = {}
	len_solution = {}
	len_target = {}

	target_categories = {}
	for target in targets:
		target_categories[target] = sbml_management.get_compounds(targets[target])

	if taxon_file:
		phylum_species, all_phylums = get_phylum(taxon_file)

	for target_categorie in target_categories:
		with open(jsons[target_categorie]) as json_data:
			dicti = json.load(json_data)
		create_stat_species(target_categories, dicti, key_species_stats_output, key_species_supdata_output, phylum_species, all_phylums)
		G = nx.Graph()
		added_node = []
		species_weight = {}
		if dicti["still_unprod"] != []:
			print("ERROR ", dicti["still_unprod"], " is unproducible")
		len_target[target_categorie] = len(dicti["newly_prod"]) + len(dicti["still_unprod"])
		len_min_sol[target_categorie] = len(dicti["bacteria"])
		len_union[target_categorie] = len(dicti["union_bacteria"])
		len_intersection[target_categorie] = len(dicti["inter_bacteria"])
		len_solution[target_categorie] = len(dicti["enum_bacteria"])
		for sol in dicti["enum_bacteria"]:
			for species_1, species_2 in combinations(dicti["enum_bacteria"][sol], 2):
				if species_1 not in added_node:
					if taxon_file:
						G.add_node(phylum_species[species_1])
					else:
						G.add_node(species_1)
					added_node.append(species_1)
				if species_2 not in added_node:
					if taxon_file:
						G.add_node(phylum_species[species_2])
					else:
						G.add_node(species_2)
					added_node.append(species_2)
				combination_species = '_'.join(sorted([species_1, species_2]))
				if combination_species not in species_weight:
					species_weight[combination_species] = 1
				else:
					species_weight[combination_species] += 1
				if taxon_file:
					G.add_edge(phylum_species[species_1], phylum_species[species_2], weight=species_weight[combination_species])
				else:
					G.add_edge(species_1, species_2, weight=species_weight[combination_species])

		with open(miscoto_stat_output,"w") as stats_output:
			statswriter = csv.writer(stats_output, delimiter='\t')
			statswriter.writerow(['categories', 'nb_target', 'size_min_sol', 'size_union', 'size_intersection', 'size_enum'])
			statswriter.writerow([target_categorie, str(len_target[target_categorie]), str(len_min_sol[target_categorie]),
								str(len_union[target_categorie]), str(len_intersection[target_categorie]), str(len_solution[target_categorie])])

		nx.write_gml(G, gml_output)

def compression(gml_input, bbl_output):
    powergrasp.compress_by_cc(gml_input, bbl_output)


def bbl_to_svg(oog_jar, bbl_input, svg_output):
    subprocess.call(['java', '-jar', oog_jar, '-inputfiles', bbl_input, '-img', '-f', svg_output])
