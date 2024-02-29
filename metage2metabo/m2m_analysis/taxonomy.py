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
import time
import sys
import logging
from ete3 import NCBITaxa

logger = logging.getLogger(__name__)


def get_taxon(taxonomy_file_path):
    """From the taxonomy file (created by extract_taxa) create a dictionary and a list linking taxon and species

    Args:
        taxonomy_file_path (str): path to the taxonomy_file

    Returns:
        taxon_named_species (dict): associate organism ID as key with taxon name as value
        all_taxa (list): list all taxa in file
    """
    taxon_named_species = {}
    all_taxa = []
    with open(taxonomy_file_path, "r") as taxonomy_file:
        taxonomy_reader = csv.reader(taxonomy_file, delimiter="\t", quotechar="|")
        next(taxonomy_reader)
        for row in taxonomy_reader:
            taxonomy = row[2].split('__')[0]
            taxon_named_species[row[0]] = row[2]
            if taxonomy not in all_taxa:
                all_taxa.append(taxonomy)

    return taxon_named_species, all_taxa


def extract_taxa(mpwt_taxon_file, taxon_output_file, tree_output_file, taxonomy_level="phylum"):
    """From NCBI taxon ID, extract taxonomy rank and create a tree file

    Args:
        mpwt_taxon_file (str): mpwt taxon file for species in sbml folder
        taxon_output_file (str): path to taxonomy output file
        tree_output_file (str): path to tree output file
        taxonomy_level (str): taxonomy level, must be: phylum, class, order, family, genus or species.
    """
    starttime = time.time()
    logger.info('######### Extract taxon information from {0}. #########'.format(mpwt_taxon_file))

    ncbi = NCBITaxa()

    # Map the taxonomy level to the taxonomy index in the ranks list.
    map_taxonomy_index = {'phylum': 0, 'class': 1, 'order': 2,
                        'family': 3, 'genus': 4, 'species': 5}
    taxonomy_index = map_taxonomy_index[taxonomy_level]

    taxon_ids = []

    taxon_count = {}
    taxonomy_file_datas = []

    with open(mpwt_taxon_file, "r") as taxon_file:
        dict_reader = csv.DictReader(taxon_file, delimiter="\t")
        headers = set(dict_reader.fieldnames)
        if set(['taxon_id', 'species']).issubset(headers):
            for line in dict_reader:
                if line['taxon_id'] is not None and line['taxon_id'] != '':
                    taxon_ids.append(line['taxon_id'])
                    lineage = ncbi.get_lineage(line['taxon_id'])
                    lineage2ranks = ncbi.get_rank(lineage)
                    names = ncbi.get_taxid_translator(lineage)
                    ranks2lineage = dict((rank, names[taxid]) for (taxid, rank) in lineage2ranks.items())
                    ranks = [ranks2lineage.get(rank, "unknown") for rank in ["phylum", "class", "order", "family", "genus", "species"]]
                    if ranks[taxonomy_index] != "unknown":
                        taxon = ranks[taxonomy_index]
                    else:
                        taxon = "unknown"

                    taxon = taxon.replace(' ', '_').replace('.', '')
                    if taxon not in taxon_count:
                        taxon_count[taxon] = 1
                    else:
                        taxon_count[taxon] += 1

                    row = ([line['species'], line['taxon_id']] + [taxon + '__' + str(taxon_count[taxon])] + ranks)
                    taxonomy_file_datas.append(row)
                else:
                    logger.warning('No taxon_id for {0} in file {1}.'.format(line['species'], mpwt_taxon_file))
        else:
            logger.critical('ERROR: No headers "taxon_id" and/or "species" in taxon file {0}.'.format(mpwt_taxon_file))
            sys.exit()

    with open(taxon_output_file, "w") as taxonomy_file:
        csvwriter = csv.writer(taxonomy_file, delimiter="\t")
        csvwriter.writerow(["organism_id", "taxid", "taxon_number", "phylum", "class", "order", "family", "genus", "species"])
        for taxonomy_file_data in taxonomy_file_datas:
            csvwriter.writerow(taxonomy_file_data)

    tree = ncbi.get_topology(taxon_ids)

    with open(tree_output_file, "w") as tree_file:
        tree_file.write(tree.get_ascii(attributes=["sci_name", "rank"]))

    logger.info(
        "--- Taxonomy runtime %.2f seconds ---\n" % (time.time() - starttime))
