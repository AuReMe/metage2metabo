import csv

from ete3 import NCBITaxa

def get_taxon(taxonomy_file_path):
    """From the taxonomy file (created by extract_taxa) create a dictionary and a list linking taxon and species

    Args:
        taxonomy_file (str): path to the taxonomy_file
    """
    taxon_named_species = {}
    all_taxons = []
    with open(taxonomy_file_path, "r") as taxonomy_file:
        taxonomy_reader = csv.reader(taxonomy_file, delimiter="\t", quotechar="|")
        next(taxonomy_reader)
        for row in taxonomy_reader:
            taxonomy = row[2].split('__')[0]
            taxon_named_species[row[0]] = row[2]
            if taxonomy not in all_taxons:
                all_taxons.append(taxonomy)

    return taxon_named_species, all_taxons

def detect_taxon_species(taxon_named_species, all_taxons):
    """From a list of species named after their taxon, return a dictionary of {taxon: [species_1, species_2]}

    Args:
        taxon_named_species (dict): {species_ID: species_named_after_taxon}
        all_taxons (list): all taxon in the dataset

    Returns:
        dict: {taxon: [species_1, species_2]}
    """
    taxon_species = {}
    for species_id in taxon_named_species:
        taxon = species_id.split('__')[0]
        if taxon not in taxon_species:
            taxon_species[taxon] = [taxon_named_species[species_id]]
        else:
            taxon_species[taxon].append(taxon_named_species[species_id])

    for taxon in all_taxons:
        if taxon not in taxon_species:
            taxon_species[taxon] = []

    return taxon_species


def extract_taxa(mpwt_taxon_file, taxon_output_file, tree_output_file, taxonomy_level="phylum"):
    """From NCBI taxon ID, extract taxonomy rank and create a tree file

    Args:
        mpwt_taxon_file (str): mpwt taxon file for species in sbml folder
        taxon_output_file (str): path to taxonomy output file
        tree_output_file (str): path to tree output file

    """
    ncbi = NCBITaxa()

    # Map the taxonomy level to the taxonomy index in the ranks list.
    map_taxonomy_index = {'phylum': 0, 'class': 1, 'order': 2,
                        'family': 3, 'genus': 4, 'species': 5}
    taxonomy_index = map_taxonomy_index[taxonomy_level]

    taxon_ids = []

    taxon_count = {}
    taxonomy_file_datas = []

    with open(mpwt_taxon_file, "r") as taxon_file:
        csvfile = csv.reader(taxon_file, delimiter="\t")
        next(csvfile)
        for line in csvfile:
            taxon_ids.append(line[1])
            lineage = ncbi.get_lineage(line[1])
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
            elif taxon == "unknown":
                taxon_count[taxon] = ""
            else:
                taxon_count[taxon] += 1

            row = ([line[0], line[1]] + [taxon + '__' + str(taxon_count[taxon])] + ranks)
            taxonomy_file_datas.append(row)

    with open(taxon_output_file, "w") as taxonomy_file:
        csvwriter = csv.writer(taxonomy_file, delimiter="\t")
        csvwriter.writerow(["organism_id", "taxid", "taxon_number", "phylum", "class", "order", "family", "genus", "species"])
        for taxonomy_file_data in taxonomy_file_datas:
            csvwriter.writerow(taxonomy_file_data)

    tree = ncbi.get_topology(taxon_ids)

    with open(tree_output_file, "w") as tree_file:
        tree_file.write(tree.get_ascii(attributes=["sci_name", "rank"]))
