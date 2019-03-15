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
import os
import mpwt


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

    args = parser.parse_args()
    inp_dir = args.genomes
    out_dir = args.output
    try:
        nb_cpu = int(args.cpu)
    except:
        print("enter a valid number of CPU")
        print(pusage)
        sys.exit(1)


    print('Hello world, I do nothing more so far ¯\_(ツ)_/¯')
    genomes_to_pgdb(inp_dir, out_dir, nb_cpu)

def genomes_to_pgdb(genomes_dir, output_dir, cpu):
    """Run Pathway Tools on each genome of the repository
    
    Args:
        genomes_dir (str): genome repository

    Returns:
        pgdb_dir (str): pgdb repository
    """
    if not os.path.isdir(genomes_dir):
        print("Genomes directory path does not exist.")
        print(pusage)
        sys.exit(1)

    pgdb_dir = output_dir + "/pgdb"
    if not utils.is_valid_dir(pgdb_dir):
        print("Impossible to access/create output directory")
        print(pusage)
        sys.exit(1)

    if not utils.check_ptools():
        print("Pathway Tools is not in the PATH, please fix it before using the program")
        print(pusage)
        sys.exit(1)

    #TODO if PGDBs are already here: prepare clean option in main and erase them if this option is selected.

    # try:
        # mpwt.multiprocess_pwt(genomes_dir, pgdb_dir, \
        #                     patho_inference=True, \
        #                     dat_creation=True, \
        #                     dat_extraction=True, \
        #                     size_reduction=False, \
        #                     number_cpu=cpu, \
        #                     patho_log=False, \
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
            print("Impossible to access/create output directory")
            print(pusage)
            sys.exit(1)

        return sbml_dir
