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
import os, os.path
import tempfile
import time
import sys
import statistics
from menetools import run_menescope
from menetools.sbml import readSBMLspecies_clyngor
from metage2metabo import utils, sbml_management
from miscoto import run_scopes, run_mincom, run_instance
from shutil import copyfile


logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)


def run_workflow(inp_dir, out_dir, nb_cpu, clean, seeds, noorphan_bool, host_mn):
    """Run the whole m2m workflow
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        nb_cpu (int): cpu number for multi-processing
        clean (bool): clean PGDB and re-run them
        seeds (str): seeds file
        noorphan_bool (bool): ignores orphan reactions if True
        host_mn (str): metabolic network file for host
    """
    # METABOLIC NETWORK RECONSTRUCTION
    sbml_dir = recon(inp_dir, out_dir, noorphan_bool, 2, nb_cpu, clean)[1]
    # INDIVIDUAL SCOPES
    union_targets_iscope = iscope(sbml_dir, seeds, out_dir)
    # COMMUNITY SCOPE
    instance_com, targets_cscope = cscope(sbml_dir, seeds, out_dir, host_mn)
    # ADDED VALUE
    newtargets = addedvalue(union_targets_iscope, targets_cscope)
    # Add these targets to the instance
    logger.info("Setting these " + str(len(newtargets)) + " as targets")
    instance_w_targets = add_targets_to_instance(
        instance_com, out_dir,
        newtargets)
    # MINCOM
    mincom(instance_w_targets, out_dir)


def recon(inp_dir, out_dir, noorphan_bool, sbml_level, nb_cpu, clean):
    """Run metabolic network reconstruction with Pathway Tools and get SBMLs
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        noorphan_bool (bool): ignores orphan reactions if True
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
                                            sbml_level, nb_cpu)
    logger.info(
        "--- Recon runtime %.2f seconds ---" % (time.time() - starttime))
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
        logger.info("--- Indiv scopes runtime %.2f seconds ---" %
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
    logger.info("--- Community scope runtime %.2f seconds ---" %
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
    logger.info("Added value of cooperation over individual metabolism: " +
                str(len(newtargets)) + " newly reachable metabolites:")
    logger.info(', '.join(newtargets))
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
    logger.info('######### Intersection of minimal communities #########')
    logger.info("# Bacteria occurring in ALL minimal community enabling to produce the target metabolites given as inputs")
    logger.info("Intersection of bacteria in minimal communities = " +
                str(len(intersection)))
    logger.info("\n".join(intersection))
    logger.info(
        "--- Mincom runtime %.2f seconds ---" % (time.time() - starttime))


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
        sys.exit(1)
    return (pgdb_dir)


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

    seed_metabolites = readSBMLspecies_clyngor(seeds, "seeds")
    logger.info("%i metabolic models considered." %(len(d_set)))
    intersection_scope = set.intersection(*list(d_set.values()))
    logger.info(str(len(intersection_scope)) + " metabolites in core reachable by all organisms (intersection)")

    union_scope = set.union(*list(d_set.values()))
    logger.info(str(len(union_scope)) + " metabolites reachable by individual organisms altogether (union), among which " + str(len(seed_metabolites)) + " seeds (growth medium)")
    len_scope = [len(d[elem]) for elem in d]
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
        for elem in target_set:
            f.write('target("' + elem + '").')
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
