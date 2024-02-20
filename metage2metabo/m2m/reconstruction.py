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
import logging
import os
import statistics
import sys
import time
import xml.etree.ElementTree as etree

from metage2metabo import utils, sbml_management

from mpwt.mpwt_workflow import multiprocess_pwt
from mpwt.utils import cleaning_input, remove_pgdbs

from multiprocessing import Pool

from padmet.classes.padmetSpec import PadmetSpec
from padmet.utils import sbmlPlugin

from menetools.sbml import get_model, get_listOfSpecies

from shutil import copyfile

logger = logging.getLogger(__name__)
logging.getLogger("mpwt").setLevel(logging.INFO)


def recon(inp_dir, out_dir, noorphan_bool, padmet_bool, sbml_level, nb_cpu, clean, use_pwt_xml):
    """Run metabolic network reconstruction with Pathway Tools and get SBMLs.
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        noorphan_bool (bool): ignores orphan reactions if True
        padmet_bool (bool): creates padmet files if True
        sbml_level (str): SBML level (2 or 3)
        nb_cpu (int): number of CPU for multiprocessing
        clean (bool): re-run metabolic reconstructions that are already available if found
        use_pwt_xml (bool): use Pathway Tools XML instead of creating them with padmet

    Returns:
        tuple: PGDB directory (str), SBML directory (str)
    """
    starttime = time.time()
    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#       Metabolic network reconstruction      #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    if use_pwt_xml and padmet_bool:
        logger.critical("-p/padmet_bool and --pwt-xml/use_pwt_xml are incompatible arguments")
        sys.exit(1)

    # Create PGDBs
    pgdb_dir = genomes_to_pgdb(inp_dir, out_dir, nb_cpu,
                                   clean, use_pwt_xml)

    if use_pwt_xml:
        sbml_dir = os.path.join(out_dir, 'sbml')
        if not os.path.exists(sbml_dir):
            os.mkdir(sbml_dir)
        for xml_file in os.listdir(pgdb_dir):
            input_xml_path = os.path.join(pgdb_dir, xml_file)
            output_xml_path = os.path.join(sbml_dir, xml_file.replace('.xml', '.sbml'))
            copyfile(input_xml_path, output_xml_path)
        padmet_folder = None

    else:
        # Create SBMLs from PGDBs
        sbml_dir = sbml_management.pgdb_to_sbml(pgdb_dir, out_dir, noorphan_bool,
                                                padmet_bool, sbml_level, nb_cpu)
        padmet_folder = os.path.join(out_dir, 'padmet')

    output_stat_file = os.path.join(out_dir, 'recon_stats.tsv')

    analyze_recon(sbml_dir, output_stat_file, padmet_folder, padmet_bool, nb_cpu)

    logger.info(
        "--- Recon runtime %.2f seconds ---\n" % (time.time() - starttime))
    return pgdb_dir, sbml_dir, padmet_folder


def genomes_to_pgdb(genomes_dir, output_dir, cpu, clean, use_pwt_xml):
    """Run Pathway Tools on each genome of the repository
    
    Args:
        genomes_dir (str): genome repository
        output_dir (str): output repository
        cpu (int): number of CPUs to use
        clean (bool): delete PGDBs in ptools-local coresponding to the input data
        use_pwt_xml (bool): use Pathway Tools XML instead of creating them with padmet

    Returns:
        pgdb_dir (str): pgdb repository
    """
    logger.info(
        "######### Running metabolic network reconstruction with Pathway Tools #########"
    )
    if not os.path.isdir(genomes_dir):
        logger.critical("Genomes directory path does not exist.")
        sys.exit(1)

    pgdb_dir = os.path.join(output_dir, 'pgdb')
    log_dir = os.path.join(output_dir,  'pgdb_log')
    ncbirc_path = os.path.join(os.path.expanduser('~'), '.ncbirc')
    log_path = os.path.join(log_dir, 'log_error.txt')

    if not utils.is_valid_dir(pgdb_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    if not utils.check_program('pathway-tools'):
        logger.critical(
            'Pathway Tools is not in the PATH, please fix it before using the program'
        )
        sys.exit(1)

    if not utils.check_program("blastp"):
        logger.critical(
            'blastp is not in the PATH, please fix it before using the program'
        )
        sys.exit(1)

    if not utils.is_valid_file(ncbirc_path):
        logger.critical(
            f'No {ncbirc_path} file, please fix it before using the program'
        )
        sys.exit(1)

    genomes_pgdbs = [genome_dir.lower() + 'cyc' for genome_dir in os.listdir(genomes_dir)]
    if clean:
        remove_pgdbs(to_delete_pgdbs=genomes_pgdbs, number_cpu=cpu)
        cleaning_input(genomes_dir, verbose=False)

    # Check whether PGDBs are already created. If yes and not --clean, pursue without running ptools again
    pgdb_dirs = [pgdb_dir.lower() + 'cyc' for pgdb_dir in os.listdir(pgdb_dir)]
    if set(pgdb_dirs) == set(genomes_pgdbs):
        logger.warning("PGDBs are already created and will be used. To overrun them, run m2m with --clean option")
        return pgdb_dir

    taxon_file = None
    if 'taxon_id.tsv' in set(next(os.walk(genomes_dir))[2]):
        taxon_file = True

    if use_pwt_xml:
        move_dat = False
        move_xml = True
    else:
        move_dat = True
        move_xml = False

    multiprocess_pwt(genomes_dir, pgdb_dir,
                        patho_inference=True,
                        patho_hole_filler=False,
                        patho_operon_predictor=False,
                        no_download_articles=False,
                        flat_creation=True,
                        dat_extraction=move_dat,
                        xml_extraction=move_xml,
                        owl_extraction=False,
                        col_extraction=False,
                        size_reduction=False,
                        number_cpu=cpu,
                        taxon_file=taxon_file,
                        patho_log=log_dir,
                        verbose=False)

    nb_genomes_dir = len([folder for folder in os.listdir(genomes_dir) if os.path.isdir(os.path.join(genomes_dir, folder))])
    if use_pwt_xml:
        nb_pgdb_dir = len([folder for folder in os.listdir(pgdb_dir) if os.path.isfile(os.path.join(pgdb_dir, folder))])
        logger.warning("Adding prefix M_ to XML form Pathway Tools.")
        for xml_file in os.listdir(pgdb_dir):
            xml_path = os.path.join(pgdb_dir, xml_file)
            update_pathway_tools_xml(xml_path, xml_path)
    else:
        nb_pgdb_dir = len([folder for folder in os.listdir(pgdb_dir) if os.path.isdir(os.path.join(pgdb_dir, folder))])

    if nb_pgdb_dir != nb_genomes_dir:
        if os.path.exists(log_path):
            logger.critical("Something went wrong running Pathway Tools. See the log file in " + log_path)
        else:
            logger.critical("Something went wrong running Pathway Tools.")
        sys.exit(1)

    return (pgdb_dir)


def create_padmet_stat(species_name, padmet_file):
    """Extract reactions/pathways/compounds/genes from a padmet file.

    Args:
        species_name (str): species names
        padmet_file (str): path to a padmet file

    Returns
        list: [species name, list of genes, list of reactions, list of reactions associated with genes, list of compounds, list of pathways]
    """
    padmetSpec = PadmetSpec(padmet_file)

    total_pwy_id = set()
    total_cpd_id = set()

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node.id for node in padmetSpec.dicOfNode.values() if node.type == "gene"]
    gene_associated_rxns = []
    genes_with_rxns = []
    rxns = []
    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        rxns.append(rxn_node.id)
        if any([rlt for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_linked_to"]):
            genes_with_rxns.extend([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_linked_to"])
            gene_associated_rxns.append(rxn_node.id)
        total_pwy_id.update(pathways_ids)

    all_pwys = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]

    genes_with_rxns = set(genes_with_rxns)

    return [species_name, genes_with_rxns, rxns, gene_associated_rxns, all_cpds, all_pwys]


def create_sbml_stat(species_name, sbml_file):
    """Extract reactions/pathways/compounds/genes from a sbml file.

    Args:
        species_name (str): species names
        sbml_file (str): path to a sbml file

    Returns
        list: [species name, list of genes, list of reactions, list of reactions associated with genes, list of compounds]
    """
    tree = etree.parse(sbml_file)
    sbml = tree.getroot()
    genes = []
    reactions = []
    gene_associated_rxns = []
    fbc_gene_associated_rxns = []
    fbc_rxn_associated_genes = []
    compounds = []
    for e in sbml:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "model":
            model_element = e
    for els in model_element:
        if 'listOfSpecies' in els.tag:
            for el in els:
                compounds.append(sbmlPlugin.convert_from_coded_id(el.get('metaid'))[0])
        if 'listOfReactions' in els.tag:
            for el in els:
                reaction_id = sbmlPlugin.convert_from_coded_id(el.get('id'))[0]
                reactions.append(reaction_id)
                for subel in el.iter():
                    if 'notes' in subel.tag:
                        for subsubel in subel.iter():
                            for subsubsubel in subsubel.iter():
                                if 'GENE_ASSOCIATION' in subsubsubel.text:
                                    for gene in sbmlPlugin.parseGeneAssoc(subsubsubel.text):
                                        if gene not in genes:
                                            genes.append(gene.replace('GENE_ASSOCIATION:', ''))
                                    if reaction_id not in gene_associated_rxns:
                                        gene_associated_rxns.append(reaction_id)
                    # Use geneProductAssociation for xml from MetaFlux.
                    elif 'geneProductAssociation' in subel.tag:
                        for subsubel in subel.iter():
                            if 'geneProductRef' in subsubel.tag:
                                gene = subsubel.get('{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct')
                                if gene:
                                    gene = gene.replace('G_', '')
                                    if gene not in fbc_rxn_associated_genes:
                                        fbc_rxn_associated_genes.append(gene)
                                    if reaction_id not in fbc_gene_associated_rxns:
                                        fbc_gene_associated_rxns.append(reaction_id)
                            else:
                                for subsubsubel in subsubel.iter():
                                    gene = subsubsubel.get('{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct')
                                    if gene:
                                        gene = gene.replace('G_', '')
                                        if gene not in fbc_rxn_associated_genes:
                                            fbc_rxn_associated_genes.append(gene)
                                        if reaction_id not in fbc_gene_associated_rxns:
                                            fbc_gene_associated_rxns.append(reaction_id)

    # For XML from MetaFlux, use genes from geneProductAssociation to get genes and reaction with genes.
    if len(genes) == 0:
        if len(fbc_rxn_associated_genes) > 0:
            genes = fbc_rxn_associated_genes

    if len(gene_associated_rxns) == 0:
        if len(fbc_gene_associated_rxns) > 0:
            gene_associated_rxns = fbc_gene_associated_rxns

    return [species_name, genes, reactions, gene_associated_rxns, compounds]


def mean_sd_data(datas):
    """Compute the mean and standard deviation from a list.

    Args:
        datas (list): list of integer/float

    Returns
        mean_data (float): mean of the list
        sd_data (flaot): standard deviation of the lsit
    """
    if len(datas) >1:
        mean_data = "{0:.2f}".format(statistics.mean(datas))
        sd_data = "(+/- {0:.2f})".format(statistics.stdev(datas))
    else:
        logger.info("No mean and standard deviation on one sample.")
        mean_data = None
        sd_data = None

    return mean_data, sd_data


def update_pathway_tools_xml(input_sbml, output_sbml):
    """ Update XML from Pathway Tools by adding a 'M_' prefix to avoid issue, when using metaboltie IDs.

    Args:
        input_sbml (str): path to xml input file
        output_sbml (str): path to xml output file
    """
    tree = etree.parse(input_sbml)
    sbml = tree.getroot()
    model = get_model(sbml)
    speciesids = [species.attrib.get("id") for species in get_listOfSpecies(model)]
    speciesids.sort(key=len, reverse=True)
    with open(input_sbml, 'r') as open_sbml:
        open_sbml_str = open_sbml.read()
        for species in speciesids:
            open_sbml_str = open_sbml_str.replace(species, 'M_'+species)

    with open(output_sbml, 'w') as open_sbml_out:
        open_sbml_out.write(open_sbml_str)


def analyze_recon(sbml_folder, output_stat_file, padmet_folder=None, padmet_bool=None, nb_cpu=1):
    """Analyze the sbml and/or the padmet files after metabolic network reconstruction.
    And write the result in a file.

    Args:
        sbml_folder (str): directory of SBML files
        output_stat_file (str): path to output stat file
        padmet_folder (str): directory of PADMET files
        padmet_bool (bool): use or not the padmet files
        nb_cpu (int): number of CPU to use
    """
    analyze_pool = Pool(processes=nb_cpu)
    if padmet_bool and padmet_folder:
        genes = {}
        reactions = {}
        gene_associated_reactions = {}
        compounds = {}
        pathways = {}

        multiprocessing_data = []

        if os.listdir(padmet_folder) == 0:
            logger.critical("No padmet in " + padmet_folder)
            sys.exit(1)

        for padmet in os.listdir(padmet_folder):
            padmet_file = os.path.join(padmet_folder, padmet)
            species_name = padmet.replace('.padmet', '')
            multiprocessing_data.append((species_name, padmet_file))
        recon_stats = analyze_pool.starmap(create_padmet_stat, multiprocessing_data)

        with open(output_stat_file, 'w') as micro_file:
            csvwriter = csv.writer(micro_file, delimiter='\t')
            csvwriter.writerow(['species', 'nb_reactions', 'nb_reactions_with_genes', 'nb_genes', 'nb_compounds', 'nb_pathways'])
            for recon_stat in recon_stats:
                species_name = recon_stat[0]
                genes[species_name] = recon_stat[1]
                reactions[species_name] = recon_stat[2]
                gene_associated_reactions[species_name] = recon_stat[3]
                compounds[species_name] = recon_stat[4]
                pathways[species_name] = recon_stat[5]
                csvwriter.writerow([species_name, len(reactions[species_name]), len(gene_associated_reactions[species_name]),
                                    len(genes[species_name]), len(compounds[species_name]), len(pathways[species_name])])
    else:
        genes = {}
        reactions = {}
        compounds = {}
        pathways = None
        gene_associated_reactions = {}

        multiprocessing_data = []

        if os.listdir(sbml_folder) == 0:
            logger.critical("No sbml in " + sbml_folder)
            sys.exit(1)

        for sbml in os.listdir(sbml_folder):
            species_name = sbml.replace('.sbml','')
            sbml_file = os.path.join(sbml_folder, sbml)
            multiprocessing_data.append((species_name, sbml_file))
        sbml_stats = analyze_pool.starmap(create_sbml_stat, multiprocessing_data)

        with open(output_stat_file, 'w') as micro_file:
            csvwriter = csv.writer(micro_file, delimiter='\t')
            csvwriter.writerow(['species', 'nb_reactions', 'nb_reactions_with_genes', 'nb_genes', 'nb_compounds'])
            for sbml_stat in sbml_stats:
                species_name = sbml_stat[0]
                genes[species_name] = set(sbml_stat[1])
                reactions[species_name] = set(sbml_stat[2])
                gene_associated_reactions[species_name] = set(sbml_stat[3])
                compounds[species_name] = set(sbml_stat[4])
                csvwriter.writerow([species_name, len(reactions[species_name]), len(gene_associated_reactions[species_name]),
                                    len(genes[species_name]), len(compounds[species_name])])

    analyze_pool.close()
    analyze_pool.join()

    logger.info("######### Stats GSMN reconstruction #########")

    if len(genes) == len(reactions) and len(genes) == len(compounds) and len(reactions) == len(compounds):
        logger.info("Number of genomes: " + str(len(genes)))

    dataset_all_reactions = set([reaction for species_name in reactions for reaction in reactions[species_name]])
    logger.info("Number of reactions in all GSMN: " + str(len(dataset_all_reactions)))

    dataset_all_compounds = set([compound for species_name in compounds for compound in compounds[species_name]])
    logger.info("Number of compounds in all GSMN: " + str(len(dataset_all_compounds)))

    species_reactions = [len(reactions[species_name]) for species_name in reactions]
    if len(species_reactions) > 1:
        mean_species_reactions, sd_species_reactions = mean_sd_data(species_reactions)
        if mean_species_reactions and sd_species_reactions:
            logger.info("Average reactions per GSMN: " + mean_species_reactions + sd_species_reactions)
    else:
        logger.info("Number of reactions in GSMN: " + str(species_reactions[0]))

    species_compounds = [len(compounds[species_name]) for species_name in compounds]
    if len(species_compounds) > 1:
        mean_species_compounds, sd_species_compounds = mean_sd_data(species_compounds)
        if mean_species_compounds and sd_species_compounds:
            logger.info("Average compounds per GSMN: " + mean_species_compounds + sd_species_compounds)
    else:
        logger.info("Number of compounds in GSMN: " + str(species_compounds[0]))

    species_genes = [len(genes[species_name]) for species_name in genes]
    if len(species_genes) > 1:
        mean_species_genes, sd_species_genes = mean_sd_data(species_genes)
        if mean_species_genes and sd_species_genes:
            logger.info("Average genes per GSMN: " + mean_species_genes + sd_species_genes)
    else:
        logger.info("Number of genes in GSMN: " + str(species_genes[0]))

    if pathways:
        species_pathways = [len(pathways[species_name]) for species_name in pathways]
        if len(species_pathways) > 1:
            mean_species_pathways, sd_species_pathways = mean_sd_data(species_pathways)
            if mean_species_pathways and sd_species_pathways:
                logger.info("Average pathways per GSMN: " + mean_species_pathways + sd_species_pathways)
        else:
            logger.info("Number of pathways in GSMN: " + str(species_pathways[0]))

    gene_reactions_assoc_percentages = []
    for species_name in reactions:
        if len(reactions[species_name]) > 0:
            gene_reactions_assoc_percentages.append(((len(gene_associated_reactions[species_name]) / len(reactions[species_name]))*100))
        else:
            gene_reactions_assoc_percentages.append(0)
            logger.info('Warning: ' + species_name + ' metabolic network contains 0 reactions.')
    if len(gene_reactions_assoc_percentages) > 1:
        mean_gene_reactions_assoc_percentages, sd_gene_reactions_assoc_percentages = mean_sd_data(gene_reactions_assoc_percentages)
        if mean_gene_reactions_assoc_percentages and sd_gene_reactions_assoc_percentages:
            logger.info('Percentage of reactions associated with genes: ' + mean_gene_reactions_assoc_percentages + sd_gene_reactions_assoc_percentages)
    else:
        logger.info('Percentage of reactions associated with genes: ' + str(gene_reactions_assoc_percentages[0]))
