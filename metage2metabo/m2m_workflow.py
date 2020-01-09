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

import json
import logging
import mpwt
import os
import tempfile
import time
import sys
import statistics
from menetools import run_menescope
from menetools.sbml import readSBMLspecies_clyngor
from metage2metabo import utils, sbml_management
from miscoto import run_scopes, run_mincom, run_instance
from shutil import copyfile
from padmet.classes.padmetSpec import PadmetSpec
import csv
import xml.etree.ElementTree as etree
from padmet.utils import sbmlPlugin
from multiprocessing import Pool


logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)


def run_workflow(inp_dir, out_dir, nb_cpu, clean, seeds, noorphan_bool, padmet_bool, host_mn):
    """Run the whole m2m workflow
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        nb_cpu (int): cpu number for multi-processing
        clean (bool): clean PGDB and re-run them
        seeds (str): seeds file
        noorphan_bool (bool): ignores orphan reactions if True
        padmet_bool (bool): creates padmet files if True
        host_mn (str): metabolic network file for host
    """
    # METABOLIC NETWORK RECONSTRUCTION
    sbml_dir = recon(inp_dir, out_dir, noorphan_bool, padmet_bool, 2, nb_cpu, clean)[1]
    # INDIVIDUAL SCOPES
    union_targets_iscope = iscope(sbml_dir, seeds, out_dir)
    # COMMUNITY SCOPE
    instance_com, targets_cscope = cscope(sbml_dir, seeds, out_dir, host_mn)
    # ADDED VALUE
    newtargets = addedvalue(union_targets_iscope, targets_cscope)
    if len(newtargets) > 0:
        # Add these targets to the instance
        logger.info("Setting these " + str(len(newtargets)) + " as targets")
        instance_w_targets = add_targets_to_instance(
            instance_com, out_dir,
            newtargets)
        # MINCOM
        mincom(instance_w_targets, out_dir)
        # remove intermediate files
        os.unlink(instance_com)
        os.unlink(instance_w_targets)
    else:
        logger.info("No newly producible compounds, hence no community selection will be computed")
        os.unlink(instance_com)


def recon(inp_dir, out_dir, noorphan_bool, padmet_bool, sbml_level, nb_cpu, clean):
    """Run metabolic network reconstruction with Pathway Tools and get SBMLs
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        noorphan_bool (bool): ignores orphan reactions if True
        padmet_bool (bool): creates padmet files if True
        nb_cpu (int): number of CPU for multiprocessing
        clean (bool): re-run metabolic reconstructions that are already available if found

    Returns:
        tuple: PGDB directory (str), SBML directory (str)
    """
    starttime = time.time()
    # Create PGDBs
    try:
        pgdb_dir = genomes_to_pgdb(inp_dir, out_dir, nb_cpu,
                                   clean)
    except:
        logger.info("Could not run Pathway Tools")
        sys.exit(1)
    # Create SBMLs from PGDBs
    sbml_dir = sbml_management.pgdb_to_sbml(pgdb_dir, out_dir, noorphan_bool,
                                            padmet_bool, sbml_level, nb_cpu)

    output_stat_file = out_dir + '/' + 'recon_stats.tsv'
    padmet_folder = out_dir + '/padmet/'
    analyze_recon(sbml_dir, output_stat_file, padmet_folder, padmet_bool, nb_cpu)

    logger.info(
        "--- Recon runtime %.2f seconds ---\n" % (time.time() - starttime))
    return pgdb_dir, sbml_dir


def iscope(sbmldir, seeds, out_dir):
    """Compute individual scopes (reachable metabolites) for SBML files in a directory
    
    Args:
        sbmldir (str): SBML files directory
        seeds (str): SBML seeds file
    
    Returns:
        set: union of reachable metabolites for all metabolic networks
    """
    # Run individual scopes of metabolic networks if any
    starttime = time.time()
    if len([
            name for name in os.listdir(sbmldir) if os.path.isfile(sbmldir + '/' + name)
            and utils.get_extension(sbmldir + '/' + name).lower() in ["xml", "sbml"]
    ]) > 1:
        scope_json = indiv_scope_run(sbmldir, seeds, out_dir)
        logger.info("Individual scopes for all metabolic networks available in " + scope_json)
        # Analyze the individual scopes results (json file)
        reachable_metabolites_union = analyze_indiv_scope(scope_json, seeds)
        logger.info("--- Indiv scopes runtime %.2f seconds ---\n" %
                    (time.time() - starttime))
        return reachable_metabolites_union
    else:
        logger.critical("The number of SBML files is <= 1 in the input directory")
        sys.exit(1)


def cscope(sbmldir, seeds, outdir, host=None):
    """Run community scope
    
    Args:
        sbmldir (str): SBML files directory
        seeds (str): SBML file for seeds
        outdir (str): output directory
        host (str, optional): Defaults to None. Host metabolic network (SBML)
    
    Returns:
        tuple: instance file (str) and community scope (set)
    """
    starttime = time.time()
    # Create instance for community analysis
    instance_com = instance_community(sbmldir, seeds, outdir, None, host)
    # Run community scope
    logger.info("Running whole-community metabolic scopes")
    community_reachable_metabolites = comm_scope_run(instance_com, outdir)
    logger.info("--- Community scope runtime %.2f seconds ---\n" %
                (time.time() - starttime))
    return instance_com, community_reachable_metabolites


def addedvalue(iscope_rm, cscope_rm):
    """Compute the added value of considering interaction with microbiota metabolism rather than individual metabolisms
    
    Args:
        iscope_rm (set): union of metabolites in all individual scopes
        cscope_rm (set): metabolites reachable by community/microbiota
    
    Returns:
        set: set of metabolites that can only be reached by a community
    """
    # Community targets = what can be produced only if cooperation occurs between species
    newtargets = cscope_rm - iscope_rm
    logger.info("\nAdded value of cooperation over individual metabolism: " +
                str(len(newtargets)) + " newly reachable metabolites: \n")
    logger.info('\n'.join(newtargets))
    logger.info("\n")
    return newtargets

def mincom(instance_w_targets, out_dir):
    """Compute minimal community selection and show analyses
    
    Args:
        instance_w_targets (str): ASP instance filepath
        out_dir (str): results directory
    """
    starttime = time.time()
    miscoto_dir = out_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)
    # Compute community selection
    logger.info("Running minimal community selection")
    all_results = compute_mincom(instance_w_targets, out_dir)

    for key in all_results:
        all_results[key] = list(all_results[key])
    with open(miscoto_dir + "/mincom.json", 'w') as dumpfile:
        json.dump(all_results, dumpfile, default=lambda x: x.__dict__)
    logger.info("Community scopes for all metabolic networks available in " +
                miscoto_dir + "/comm_scopes.json")
    # Give one solution
    one_sol_bact = []
    for bact in all_results['bacteria']:
        one_sol_bact.append(bact)
    logger.info('######### One minimal community #########')
    logger.info("# One minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Minimal number of bacteria in communities = " +
                str(len(one_sol_bact)))
    logger.info("\n".join(one_sol_bact))
    # Give union of solutions
    union = all_results['union_bacteria']
    logger.info('######### Keystone species: Union of minimal communities #########')
    logger.info("# Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Keystone species = " +
                str(len(union)))
    logger.info("\n".join(union))
    # Give intersection of solutions
    intersection = all_results['inter_bacteria']
    logger.info('######### Essential symbionts: Intersection of minimal communities #########')
    logger.info("# Bacteria occurring in ALL minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Essential symbionts = " +
                str(len(intersection)))
    logger.info("\n".join(intersection))
    # Give keystones, essential and alternative symbionts
    alternative_symbionts = list(set(union) - set(intersection))
    logger.info('######### Alternative symbionts: Difference between Union and Intersection #########')
    logger.info("# Bacteria occurring in at least one minimal community but not all minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Alternative symbionts = " +
                str(len(alternative_symbionts)))
    logger.info("\n".join(alternative_symbionts))
    logger.info(
        "--- Mincom runtime %.2f seconds ---\n" % (time.time() - starttime))


def genomes_to_pgdb(genomes_dir, output_dir, cpu, clean):
    """Run Pathway Tools on each genome of the repository
    
    Args:
        genomes_dir (str): genome repository

    Returns:
        pgdb_dir (str): pgdb repository
    """
    logger.info(
        "######### Running metabolic network reconstruction with Pathway Tools #########"
    )
    if not os.path.isdir(genomes_dir):
        logger.critical("Genomes directory path does not exist.")
        sys.exit(1)

    pgdb_dir = output_dir + "/pgdb"
    if not utils.is_valid_dir(pgdb_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    if not utils.check_program("pathway-tools"):
        logger.critical(
            "Pathway Tools is not in the PATH, please fix it before using the program"
        )
        sys.exit(1)

    if not utils.check_program("blastp"):
        logger.critical(
            "blastp is not in the PATH, please fix it before using the program"
        )
        sys.exit(1)

    if not utils.is_valid_file(os.path.expanduser("~") + "/.ncbirc"):
        logger.critical(
            "No ~/.ncbirc file, please fix it before using the program"
        )
        sys.exit(1)

    genomes_pgdbs = [genome_dir.lower() + 'cyc' for genome_dir in os.listdir(genomes_dir)]
    already_here_pgdbs = mpwt.list_pgdb()
    if clean and set(genomes_pgdbs).issubset(set(already_here_pgdbs)):
        mpwt.cleaning(number_cpu=cpu, verbose=True)

    # Check whether PGDBs are already created. If yes and not --clean, pursue without running ptools again
    pgdb_dirs = [pgdb_dir.lower() + 'cyc' for pgdb_dir in os.listdir(pgdb_dir)]
    if set(pgdb_dirs) ==set(genomes_pgdbs):
        logger.warning("PGDBs are already created and will be used. To overrun them, run m2m with --clean option")
        return pgdb_dir

    try:
        mpwt.multiprocess_pwt(genomes_dir, pgdb_dir,
                            patho_inference=True,
                            patho_hole_filler=False,
                            patho_operon_predictor=False,
                            no_download_articles=False,
                            dat_creation=True,
                            dat_extraction=True,
                            size_reduction=False,
                            number_cpu=cpu,
                            patho_log=False,
                            verbose=True)
    except:
        logger.critical("Something went wrong running Pathway Tools")
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
    rxns = []
    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        rxns.append(rxn_node.id)
        if any([rlt for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_linked_to"]):
            gene_associated_rxns.append(rxn_node.id)
        total_pwy_id.update(pathways_ids)

    all_pwys = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]

    return [species_name, all_genes, rxns, gene_associated_rxns, all_cpds, all_pwys]


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
                for subel in el.getchildren():
                    if 'notes' in subel.tag:
                        for subsubel in subel.getchildren():
                            for subsubsubel in subsubel.getchildren():
                                if 'GENE_ASSOCIATION' in subsubsubel.text:
                                    for gene in sbmlPlugin.parseGeneAssoc(subsubsubel.text):
                                        genes.append(gene.replace('GENE_ASSOCIATION:', ''))
                                    if reaction_id not in gene_associated_rxns:
                                        gene_associated_rxns.append(reaction_id)

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

def analyze_recon(sbml_folder, output_stat_file, padmet_folder, padmet_bool=None, nb_cpu=1):
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
    if padmet_bool:
        genes = {}
        reactions = {}
        gene_associated_reactions = {}
        compounds = {}
        pathways = {}

        multiprocessing_data = []
        for padmet in os.listdir(padmet_folder):
            padmet_file = padmet_folder + '/' + padmet
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
        for sbml in os.listdir(sbml_folder):
            species_name = sbml.replace('.sbml','')
            sbml_file = sbml_folder + '/' + sbml
            multiprocessing_data.append((species_name, sbml_file))
        sbml_stats = analyze_pool.starmap(create_sbml_stat, multiprocessing_data)

        with open(output_stat_file, 'w') as micro_file:
            csvwriter = csv.writer(micro_file, delimiter='\t')
            csvwriter.writerow(['species', 'nb_reactions', 'nb_genes', 'nb_compounds'])
            for sbml_stat in sbml_stats:
                species_name = sbml_stat[0]
                genes[species_name] = set(sbml_stat[1])
                reactions[species_name] = set(sbml_stat[2])
                gene_associated_reactions[species_name] = set(sbml_stat[3])
                compounds[species_name] = set(sbml_stat[4])
                csvwriter.writerow([species_name, len(reactions[species_name]), len(genes[species_name]), len(compounds[species_name])])

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
        gene_reactions_assoc_percentages.append(((len(gene_associated_reactions[species_name]) / len(reactions[species_name]))*100))
    if len(gene_reactions_assoc_percentages) > 1:
        mean_gene_reactions_assoc_percentages, sd_gene_reactions_assoc_percentages = mean_sd_data(gene_reactions_assoc_percentages)
        if mean_gene_reactions_assoc_percentages and sd_gene_reactions_assoc_percentages:
            logger.info('Percentage of reactions associated with genes: ' + mean_gene_reactions_assoc_percentages + sd_gene_reactions_assoc_percentages)
    else:
        logger.info('Percentage of reactions associated with genes: ' + str(gene_reactions_assoc_percentages[0]))


def indiv_scope_run(sbml_dir, seeds, output_dir):
    """Run Menetools and analyse individual metabolic capabilities
    
    Args:
        sbml_dir (str): directory of SBML files
        output_dir (str): directory for results
    
    Returns:
        str: output file for Menetools analysis
    """
    logger.info("######### Running individual metabolic scopes #########")
    menetools_dir = output_dir + "/indiv_scopes"
    if not utils.is_valid_dir(menetools_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    all_files = [
        f for f in os.listdir(sbml_dir)
        if os.path.isfile(os.path.join(sbml_dir, f)) and utils.get_extension(
            os.path.join(sbml_dir, f)).lower() in ["xml", "sbml"]
    ]
    all_scopes = {}
    for f in all_files:
        bname = utils.get_basename(f)
        try:
            all_scopes[bname] = run_menescope(
                draft_sbml=os.path.join(sbml_dir, f), seeds_sbml=seeds)
        except:
            logger.critical("Something went wrong running Menetools")

    with open(menetools_dir + "/indiv_scopes.json", 'w') as dumpfile:
        json.dump(all_scopes, dumpfile)
    return menetools_dir + "/indiv_scopes.json"


def analyze_indiv_scope(jsonfile, seeds):
    """Analyze the output of Menescope, stored in a json
    
    Args:
        jsonfile (str): output of menescope
    """
    with open(jsonfile) as json_data:
        d = json.load(json_data)
        if not d:
            logger.critical("Json file is empty. Individual scopes calculation failed. Please fill an issue on Github")
            sys.exit(1)
    d_set = {}

    for elem in d:
        d_set[elem] = set(d[elem])

    try:
        seed_metabolites = readSBMLspecies_clyngor(seeds, "seeds")
    except FileNotFoundError:
        logger.critical("File not found: "+seeds)
        sys.exit(1)
    except etree.ParseError:
        logger.critical("Invalid syntax in SBML file: "+seeds)
        sys.exit(1)

    logger.info("%i metabolic models considered." %(len(d_set)))
    intersection_scope = set.intersection(*list(d_set.values()))
    logger.info('\n' + str(len(intersection_scope)) + " metabolites in core reachable by all organisms (intersection) \n")
    logger.info("\n".join(intersection_scope))

    union_scope = set.union(*list(d_set.values()))
    logger.info('\n' + str(len(union_scope)) + " metabolites reachable by individual organisms altogether (union), among which " + str(len(seed_metabolites)) + " seeds (growth medium) \n")
    logger.info("\n".join(union_scope))
    len_scope = [len(d[elem]) for elem in d]
    logger.info("\nintersection of scope " + str(len(intersection_scope)))
    logger.info("union of scope " + str(len(union_scope)))
    logger.info("max metabolites in scope " + str(max(len_scope)))
    logger.info("min metabolites in scope " + str(min(len_scope)))
    logger.info("average number of metabolites in scope %.2f (+/- %.2f)" %
                (statistics.mean(len_scope), statistics.stdev(len_scope)))
    return union_scope


def instance_community(sbml_dir, seeds, output_dir, targets = None, host_mn=None):
    """Create ASP instance for community analysis
    
    Args:
        sbml_dir (str): directory of symbionts SBML files
        seeds (str): seeds SBML file
        output_dir (str): directory for results

    Returns:
        str: instance filepath
    """
    logger.info(
            "######### Creating metabolic instance for the whole community #########"
        )
    miscoto_dir = output_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    # create tempfile
    fd, outputfile = tempfile.mkstemp(suffix='.lp', prefix='miscoto_', dir=miscoto_dir)

    instance_filepath = run_instance(
        bacteria_dir=sbml_dir,
        seeds_file=seeds,
        host_file=host_mn,
        targets_file=targets,
        output=outputfile)

    logger.info("Created instance in " + instance_filepath)
    return instance_filepath


def comm_scope_run(instance, output_dir):
    """Run Miscoto_scope and analyse individual metabolic capabilities
    
    Args:
        sbml_dir (str): directory of SBML files
        seeds (str): seeds SBML file
        output_dir (str): directory for results
    
    Returns:
        set: microbiota scope
    """
    miscoto_dir = output_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)
    microbiota_scope = run_scopes(instance)
    with open(miscoto_dir + "/comm_scopes.json", 'w') as dumpfile:
        json.dump(microbiota_scope, dumpfile)
    logger.info("Community scopes for all metabolic networks available in " +
                miscoto_dir + "/comm_scopes.json")
    return set(microbiota_scope['com_scope'])


def add_targets_to_instance(instancefile, output_dir, target_set):
    """Add targets to the ASP community instance file
    
    Args:
        instancefile (str): instance filepath
        target_set (set): targets to be added
    
    Returns:
        str: new instance filepath
    """
    new_instance_file = output_dir + "/community_analysis/" + utils.get_basename(instancefile) + '__tgts.lp'
    copyfile(instancefile, new_instance_file)

    with open(new_instance_file, 'a') as f:
        f.write('\n')
        for elem in target_set:
            f.write('target("' + elem + '").\n')
    return new_instance_file


def compute_mincom(instancefile, output_dir):
    """Run minimal community selection and analysis
    
    Args:
        instancefile (str): filepath to instance file
        output_dir (str): directory with results
    
    Returns:
        dict: results of miscoto_mincom analysis
    """
    miscoto_dir = output_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    results_dic = run_mincom(option="soup",
                            lp_instance_file=instancefile,
                            optsol=True,
                            union=True,
                            intersection=True)
    return results_dic
