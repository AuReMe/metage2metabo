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
import time
from metage2metabo import utils, padmet2sbml
import os, os.path
import mpwt
from menetools import run_menescope
import json
import logging

logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)

###############################################################################
#
message = """
Description here
"""

pusage = """
Usage here
"""

requires = """
Requirements here"
"""
#
###############################################################################


def run_workflow():
    """description
    """
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

    args = parser.parse_args()
    inp_dir = args.genomes
    out_dir = args.output
    seeds = args.seeds
    try:
        nb_cpu = int(args.cpu)
    except:
        logger.critical("enter a valid number of CPU")
        logger.info(pusage)
        sys.exit(1)

    plop = indiv_scope_run(
        "/Users/cfrioux/wd/scripts/metage2metabo/toy_res/sbml",
        "/Users/cfrioux/wd/ecosystems/bft_and_co/Seafile/sync_stage_enora/data/alga/seeds_modified_for_ions.xml",
        "/Users/cfrioux/wd/scripts/metage2metabo/toy_res/")
    analyze_indiv_scope(plop)
    # print('Hello world, I do nothing more so far ¯\_(ツ)_/¯')
    # genomes_to_pgdb(inp_dir, out_dir, nb_cpu)

def genomes_to_pgdb(genomes_dir, output_dir, cpu):
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

    #TODO if PGDBs are already here: prepare clean option in main and erase them if this option is selected.

    # try:
    # mpwt.multiprocess_pwt(genomes_dir, pgdb_dir,
    #                     patho_inference=True,
    #                     dat_creation=True,
    #                     dat_extraction=True,
    #                     size_reduction=False,
    #                     number_cpu=cpu,
    #                     patho_log=False,
    #                     verbose=True)
    # except:
    #     print("Oops, something went wrong running Pathway Tools")

    return (pgdb_dir)

def pgdb_to_sbml(pgdb_dir, output_dir):
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

    return sbml_dir

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
        # try:
        all_scopes[bname] = run_menescope(
            draft_sbml=os.path.join(sbml_dir, f), seeds_sbml=seeds)
        # except:
        #     print("Oops, something went wrong running Menetools")

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
    logger.info(str(len(intersection_scope)) + " metabolites in intersection")

    union_scope = set.union(*list(d_set.values()))
    len_scope = [len(d[elem]) for elem in d]
    logger.info("max metabolites in scope " + str(max(len_scope)))
    logger.info("min metabolites in scope " + str(min(len_scope)))
    return


if __name__ == "__main__":
    run_workflow()