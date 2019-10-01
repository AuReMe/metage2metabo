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

recipes = '''{
	"clingo multithreading": 4
}'''
powergrasp.recipe.Recipe.from_lines(recipes)

from ete3 import NCBITaxa
from itertools import combinations
from metage2metabo import utils, sbml_management


def run_analysis(m2m_results_folder, target_file, seed_file, output_dir, taxon_file):
	sbml_folder = m2m_results_folder + '/sbml'
	miscoto_analysis = output_dir + '/miscoto_analysis/'

	results = miscoto.run_mincom(option='soup', bacteria_dir=sbml_folder,
						targets_file=target_file, seeds_file=seed_file, host_file=None,
                		intersection=True, enumeration=True, union=True, optsol=True)

	miscoto_json = miscoto_analysis + 'mincom.json'
	miscoto.utils.results_to_json(results, miscoto_json)

	if not utils.is_valid_dir(miscoto_analysis):
		logger.critical("Impossible to access/create output directory")
		sys.exit(1)

	miscoto_stat_output = miscoto_analysis + 'miscoto_stats.txt'
	gml_output = miscoto_analysis + 'graph.gml'
	bbl_output = miscoto_analysis + 'graph.bbl'

	create_gml(sbml_folder, target_file, miscoto_json, miscoto_stat_output, gml_output, taxon_file)
	compression(gml_output, bbl_output, recipes)


def extract_taxa(mpwt_taxon_file, taxon_output_file, tree_output_file):
	ncbi = NCBITaxa()

	taxon_ids = []

	phylum_count = {}
	with open(taxon_output_file, 'w') as phylum_file:
		csvwriter = csv.writer(phylum_file, delimiter='\t')
		csvwriter.writerow(['species', 'taxid', 'ink_name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
		with open(mpwt_taxon_file, 'r') as taxon_file:
			csvfile = csv.reader(taxon_file, delimiter='\t')
			for line in csvfile:
				if 'taxon' not in line[1]:
					taxon_ids.append(line[1])
					lineage = ncbi.get_lineage(line[1])
					lineage2ranks = ncbi.get_rank(lineage)
					names = ncbi.get_taxid_translator(lineage)
					ranks2lineage = dict((rank, names[taxid]) for (taxid, rank) in lineage2ranks.items())
					ranks = [ranks2lineage.get(rank, '<not present>') for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]
					phylum = ranks[1][:4]
					if phylum not in phylum_count :
						phylum_count[phylum] = 1
					else:
						phylum_count[phylum] += 1
					row = [line[0], line[1]] + [phylum + str(phylum_count[phylum])] + ranks
					csvwriter.writerow(row)

	tree = ncbi.get_topology(taxon_ids)

	with open(tree_output_file, 'w') as tree_file:
		tree_file.write(tree.get_ascii(attributes=["sci_name", "rank"]))


def create_gml(sbml_folder, target_file, json_miscoto_file, miscoto_stat_output, gml_output, taxon_file):
	len_min_sol = {}
	len_union = {}
	len_intersection = {}
	len_solution = {}
	len_target = {}

	all_species = []
	for species in os.listdir(sbml_folder):
		species = species.replace('.sbml', '')
		all_species.append(species)

	species_categories = {}
	species_categories['target'] = sbml_management.get_compounds(target_file)

	if taxon_file:
		phylum_species = {}
		with open(taxon_file, 'r') as phylum_file:
			phylum_reader = csv.reader(phylum_file, delimiter='\t', quotechar='|')
			for row in phylum_reader:
				phylum_species[row[0]] = row[2]

	for categories in  species_categories:
		with open(json_miscoto_file) as json_data:
			dicti = json.load(json_data)
		G = nx.Graph()
		added_node = []
		species_weight = {}
		if dicti["still_unprod"] != []:
			print("ERROR ", dicti["still_unprod"], " is unproducible")
		len_target[categories] = len(dicti["newly_prod"]) + len(dicti["still_unprod"])
		len_min_sol[categories] = len(dicti["bacteria"])
		len_union[categories] = len(dicti["union_bacteria"])
		len_intersection[categories] = len(dicti["inter_bacteria"])
		len_solution[categories] = len(dicti["enum_bacteria"])
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
			statswriter.writerow([categories, str(len_target[categories]), str(len_min_sol[categories]), str(len_union[categories]), str(len_intersection[categories]), str(len_solution[categories])])

		nx.write_gml(G, gml_output)

def compression(gml_input, bbl_output, recipes):
    powergrasp.compress_by_cc(gml_input, bbl_output)

def main_analysis():
	"""Run programm
	"""
	start_time = time.time()
	parser = argparse.ArgumentParser(
		"m2m_analysis",
		description=" Run graph creation and compression on enumeration of miscoto solution."
	)

	parser.add_argument(
		"-n",
		"--networksdir",
		metavar="NETWORKS_DIR",
		help="metabolic networks directory",
		required=True)

	parser.add_argument(
		"-o",
		"--out",
		dest="out",
		required=True,
		help="output directory path",
		metavar="OUPUT_DIR")

	parser.add_argument(
		"-s",
		"--seeds",
		help="seeds (growth medium) for metabolic analysis",
		required=True)

	parser.add_argument(
		"-t",
		"--targets",
		help="targets for metabolic analysis",
		required=True)

	parser.add_argument(
		"-j",
		"--json",
		help="json from miscoto",
		required=True)

	parser.add_argument(
		"-m",
		"--modelhost",
		help="host metabolic model for community analysis",
		required=False,
		default=None)

	parser.add_argument(
		"-q",
		"--quiet",
		dest="quiet",
		help="quiet mode",
		required=False,
		action="store_true",
		default=None,
	)

	parser.add_argument(
		"-c",
		"--cpu",
		help="cpu number for multi-process",
		required=False,
		type=int,
		default=1)

	parser.add_argument(
		"--taxon",
		metavar="TAXON_FILE",
		help="mpwt taxon file tsv", required=False)

    # If no argument print the help.
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	run_analysis(args.networksdir, args.targets, args.seeds, args.out, args.taxon)