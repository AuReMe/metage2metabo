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
import json
import logging
import mpwt
import os, os.path
import pyasp
import tempfile
import time
import sys

from menetools import run_menescope
from metage2metabo import utils
from miscoto import run_scopes, run_mincom, run_instance
from multiprocessing import Pool
from padmet_utils.scripts.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet_utils.scripts.connection.sbmlGenerator import padmet_to_sbml
from padmet_utils.scripts.connection.sbml_to_sbml import from_sbml_to_sbml
from padmet_utils.scripts.connection.sbml_to_padmet import from_sbml_to_padmet
from shutil import copyfile

logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)

###############################################################################
#
message = """
Description here TODO
"""

pusage = """
Usage here TODO
"""

requires = """
Requirements here TODO
"""
#
###############################################################################


def run_workflow():
    """description
    """
    # Check pyasp binaries.
    pyasp_bin_path = pyasp.__path__[0] + '/bin/'
    bin_check = []
    for asp_bin in ['clasp', 'gringo3', 'gringo4']:
        bin_check.append(utils.is_valid_path(pyasp_bin_path + asp_bin))

    if not all(bin_check):
        logger.critical("Error with pyasp installation, retry to install it:")
        logger.critical("pip install pyasp==1.4.3 --no-cache-dir --force-reinstall")
        logger.info(pusage)
        sys.exit(1)

    parser = argparse.ArgumentParser(description=message, usage=pusage, epilog=requires)
    parser.add_argument("-g",
                        "--genomes",
                        help="annotated genomes directory",
                        required=True)
    parser.add_argument("-o",
                        "--output",
                        help="output directory",
                        required=True)
    parser.add_argument("-c",
                        "--cpu",
                        help="cpu number for multi-process Pathway Tools",
                        required=False,
                        default=1)
    parser.add_argument("-s",
                        "--seeds",
                        help="seeds (growth medium) for metabolic analysis",
                        required=True)
    parser.add_argument("-m",
                        "--modelhost",
                        help="host metabolic model for community analysis",
                        required=False,
                        default=None)
    parser.add_argument("--clean",
                        help="clean PGDBs if already present",
                        required=False,
                        action="store_true")

    args = parser.parse_args()
    inp_dir = args.genomes
    out_dir = args.output
    seeds = args.seeds
    host = args.modelhost
    clean = args.clean

    try:
        nb_cpu = int(args.cpu)
    except:
        logger.critical("enter a valid number of CPU")
        logger.info(pusage)
        sys.exit(1)

    # METABOLIC NETWORK RECONSTRUCTION
    # Create PGDBs
    logger.info("######### Running metabolic network reconstruction with Pathway Tools #########")
    try:
        pgdb_dir = genomes_to_pgdb(inp_dir, out_dir, nb_cpu, clean)
    except:
        logger.info("Could not run Pathway Tools")
        sys.exit(1)
    # Create SBMLs from PGDBs
    logger.info("######### Creating SBML files #########")
    sbml_dir = pgdb_to_sbml(pgdb_dir, out_dir, nb_cpu)
    # ANALYSIS
    # Run individual scopes of metabolic networks if any
    if len([
            name for name in os.listdir(sbml_dir) if os.path.isfile(sbml_dir + '/' + name)
            and utils.get_extension(sbml_dir + '/' + name).lower() in ["xml", "sbml"]
    ]) > 1:
        logger.info("######### Running individual metabolic scopes #########")
        scope_json = indiv_scope_run(sbml_dir, seeds, out_dir)
        # Analyze the individual scopes results (json file)
        uniontargets = analyze_indiv_scope(scope_json)
        # Create instance for community analysis
        logger.info(
            "######### Creating metabolic instance for the whole community #########"
        )
        instance_com = instance_community(sbml_dir, seeds, out_dir)
        # Run community scope
        logger.info("Running whole-community metabolic scopes")
        microbiotascope = comm_scope_run(instance_com, out_dir)
        # Community targets = what can be produced only if cooperation occurs between species
        newtargets = set(microbiotascope) - uniontargets
        logger.info(str(len(newtargets)) + " metabolites can only be produced through metabolic cooperation")
        logger.info(newtargets)
        logger.info("Setting these " + str(len(newtargets)) + " as targets")
        # Add these targets to the instance
        instance_w_targets = add_targets_to_instance(
            instance_com, out_dir,
            newtargets)
        # Compute community selection
        logger.info("Running minimal community selection")
        all_results = mincom(instance_w_targets, out_dir)
        # Give one solution
        onesol = all_results['one_model']
        one_sol_bact = []
        for a in onesol:
            if a.pred() == 'chosen_bacteria':
                one_sol_bact.append(a.arg(0).rstrip('"').lstrip('"'))
        logger.info('######### One minimal community #########')
        logger.info("# One minimal community enabling to produce the target metabolites given as inputs")
        logger.info("Minimal number of bacteria in communities = " +
                    str(len(one_sol_bact)))
        logger.info("\n".join(one_sol_bact))
        # Give union of solutions
        union = all_results['union_bacteria']
        logger.info('######### Union of minimal communities #########')
        logger.info("# Bacteria occurring in at least one minimal community enabling to produce the target metabolites given as inputs")
        logger.info("Union of bacteria in minimal communities = " +
                    str(len(union)))
        logger.info("\n".join(union))
        # Give intersection of solutions
        intersection = all_results['inter_bacteria']
        logger.info('######### Union of minimal communities #########')
        logger.info("# Bacteria occurring in ALL minimal community enabling to produce the target metabolites given as inputs")
        logger.info("Intersection of bacteria in minimal communities = " +
                    str(len(intersection)))
        logger.info("\n".join(intersection))


def genomes_to_pgdb(genomes_dir, output_dir, cpu, clean):
    """Run Pathway Tools on each genome of the repository
    
    Args:
        genomes_dir (str): genome repository

    Returns:
        pgdb_dir (str): pgdb repository
    """
    if not os.path.isdir(genomes_dir):
        logger.critical("Genomes directory path does not exist.")
        logger.info(pusage)
        sys.exit(1)

    pgdb_dir = output_dir + "/pgdb"
    if not utils.is_valid_dir(pgdb_dir):
        logger.critical("Impossible to access/create output directory")
        logger.info(pusage)
        sys.exit(1)

    if not utils.check_ptools():
        logger.critical(
            "Pathway Tools is not in the PATH, please fix it before using the program"
        )
        logger.info(pusage)
        sys.exit(1)

    #TODO test whether blast is in the PATH

    #TODO test is ncbirc is in home folder

    #TODO test whether not everything is run again in case some PGDBs already are in PDGB_dir and/or Ptools local

    genomes_pgdbs = [genome_dir.lower() + 'cyc' for genome_dir in os.listdir(genomes_dir)]
    already_here_pgdbs = mpwt.list_pgdb()
    if clean and set(genomes_pgdbs).issubset(set(already_here_pgdbs)):
        mpwt.cleaning(cpu)

    # Check whether PGDBs are already created. If yes and not --clean, pursue without running ptools again
    pgdb_dirs = [pgdb_dir.lower() + 'cyc' for pgdb_dir in os.listdir(pgdb_dir)]
    if set(pgdb_dirs) ==set(genomes_pgdbs):
        logger.warning("PGDBs are already created and will be used. To overrun them, run m2m with --clean option")
        return pgdb_dir

    try:
        mpwt.multiprocess_pwt(genomes_dir, pgdb_dir,
                            patho_inference=True,
                            dat_creation=True,
                            dat_extraction=True,
                            size_reduction=False,
                            number_cpu=cpu,
                            patho_log=False,
                            verbose=True)
    except:
        logger.critical("Something went wrong running Pathway Tools")

    return (pgdb_dir)

def run_pgdb_to_sbml(species_multiprocess_data):
    """Function used in multiprocess to turn PGDBs into SBML2.

    Args:
        species_multiprocess_data (list): pathname to species pgdb dir, pathname to species sbml file

    Returns:
        sbml_check (bool): Check if sbml file exists
    """
    species_pgdb_dir = species_multiprocess_data[0]
    species_sbml_file = species_multiprocess_data[1]
    padmet = from_pgdb_to_padmet(species_pgdb_dir, 'metacyc', '22.5', arg_verbose=False, arg_with_genes=True, arg_source=None, enhanced_db=None, padmetRef_file=None)

    padmet_to_sbml(padmet, species_sbml_file, sbml_lvl=2, verbose=False)

    sbml_check = utils.is_valid_path(species_sbml_file)

    return sbml_check

def pgdb_to_sbml(pgdb_dir, output_dir, cpu):
    """Turn Pathway Tools PGDBs into SBML2 files using Padmet
    
    Args:
        pgdb_dir (str): PGDB directory

    Returns:
        sbml_dir (str): SBML directory
    """
    sbml_dir = output_dir + "/sbml"
    if not utils.is_valid_dir(sbml_dir):
        logger.critical("Impossible to access/create output directory")
        logger.info(pusage)
        sys.exit(1)

    pgdb_to_sbml_pool = Pool(processes=cpu)

    multiprocess_data = []
    for species in os.listdir(pgdb_dir):
        multiprocess_data.append([pgdb_dir + '/' + species, sbml_dir + '/' + species + '.sbml'])

    sbml_checks = pgdb_to_sbml_pool.map(run_pgdb_to_sbml, multiprocess_data)

    pgdb_to_sbml_pool.close()
    pgdb_to_sbml_pool.join()

    if all(sbml_checks):
        return sbml_dir

    else:
        logger.critical("Error during padmet/sbml creation.")
        logger.info(pusage)
        sys.exit(1)

def sbml3_to_sbml2(sbml_file, sbml_output_file, cpu, version):
    """Turn sbml3 to sbml2.

    Args:
        sbml_file (string): pathname to species sbml file
        sbml_output_file (string): pathname to output sbml file
        db (string): name of the database ('metacyc')
        version (string): version of tthe database

    Returns:
        sbml_check (bool): Check if sbml file exists
    """
    from_sbml_to_sbml(sbml_file, sbml_output_file, version, cpu)

    sbml_check = utils.is_valid_path(sbml_output_file)

    return sbml_check

def sbml2_to_sbml3(sbml_file, sbml_output_file, db, version):
    """Turn sbml2 to sbml3.

    Args:
        sbml_file (string): pathname to species sbml file
        sbml_output_file (string): pathname to output sbml file
        db (string): name of the database ('metacyc')
        version (string): version of tthe database

    Returns:
        sbml_check (bool): Check if sbml file exists
    """
    padmet = from_sbml_to_padmet(sbml=sbml_file, db=db, version=version, padmetSpec_file=None, source_tool=None, source_category=None, source_id=None, padmetRef_file=None, mapping=None, verbose=None)
    padmet_to_sbml(padmet, sbml_output_file, sbml_lvl=3, verbose=False)

    sbml_check = utils.is_valid_path(sbml_output_file)

    return sbml_check

def indiv_scope_run(sbml_dir, seeds, output_dir):
    """Run Menetools and analyse individual metabolic capabilities
    
    Args:
        sbml_dir (str): directory of SBML files
        output_dir (str): directory for results
    
    Returns:
        str: output file for Menetools analysis
    """
    menetools_dir = output_dir + "/indiv_scopes"
    if not utils.is_valid_dir(menetools_dir):
        logger.critical("Impossible to access/create output directory")
        logger.info(pusage)
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
            #TODO catch OSError and look for ASP binaries  pip install pyasp==1.4.3 --no-cache-dir --force-reinstall
            logger.critical("Something went wrong running Menetools")

    with open(menetools_dir + "/indiv_scopes.json", 'w') as dumpfile:
        json.dump(all_scopes, dumpfile)

    return menetools_dir + "/indiv_scopes.json"


def analyze_indiv_scope(jsonfile):
    """Analyze the output of Menescope, stored in a json
    
    Args:
        jsonfile (str): output of menescope
    """
    with open(jsonfile) as json_data:
        d = json.load(json_data)
    d_set = {}

    for elem in d:
        d_set[elem] = set(d[elem])

    intersection_scope = set.intersection(*list(d_set.values()))
    logger.info(str(len(intersection_scope)) + " metabolites in core reachable by all organisms (intersection)")

    union_scope = set.union(*list(d_set.values()))
    logger.info(str(len(union_scope)) + " metabolites reachable by individual organisms altogether (union)")
    len_scope = [len(d[elem]) for elem in d]
    logger.info("max metabolites in scope " + str(max(len_scope)))
    logger.info("min metabolites in scope " + str(min(len_scope)))

    return union_scope

def instance_community(sbml_dir, seeds, output_dir, host_mn=None):
    """create ASP instance for community analysis
    
    Args:
        sbml_dir (str): directory of symbionts SBML files
        seeds (str): seeds SBML file
        output_dir (str): directory for results

    Returns:
        str: instance filepath
    """
    miscoto_dir = output_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        logger.info(pusage)
        sys.exit(1)

    # create tempfile
    fd, outputfile = tempfile.mkstemp(suffix='.lp', prefix='miscoto_', dir=miscoto_dir)

    instance_filepath = run_instance(
        bacteria_dir=sbml_dir,
        seeds_file=seeds,
        host_file=host_mn,
        targets_file=None,
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
        lst: microbiota scope
    """
    miscoto_dir = output_dir + "/community_analysis"
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        logger.info(pusage)
        sys.exit(1)

    microbiota_scope = run_scopes(instance)
    return microbiota_scope['com_scope']

def add_targets_to_instance(instancefile, output_dir, target_set):
    """Add targets to the ASP community instance file
    
    Args:
        instancefile (str): instance filepath
        target_set (set): targets to be added
    
    Returns:
        str: new instance filepath
    """
    new_instance_file = output_dir + "/community_analysis" + utils.get_basename(instancefile) + '__tgts.lp'
    copyfile(instancefile, new_instance_file)

    with open(new_instance_file, 'a') as f:
        for elem in target_set:
            f.write('target("' + elem + '").')

    return new_instance_file

def mincom(instancefile, output_dir):
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
        logger.info(pusage)
        sys.exit(1)

    results_dic = run_mincom(option="soup",
                            lp_instance_file=instancefile,
                            optsol=True,
                            union=True,
                            intersection=True)
    return results_dic

if __name__ == "__main__":
    run_workflow()