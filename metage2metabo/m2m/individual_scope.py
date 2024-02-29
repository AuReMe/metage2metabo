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
import json
import logging
import os
import statistics
import sys
import time
import traceback
import xml.etree.ElementTree as etree

from menetools import run_menescope
from menetools.sbml import readSBMLspecies_clyngor

from metage2metabo import utils

from multiprocessing import Pool

logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)


def iscope(sbmldir, seeds, out_dir, cpu_number=1):
    """Compute individual scopes (reachable metabolites) for SBML files in a directory.
    
    Args:
        sbmldir (str): SBML files directory
        seeds (str): SBML seeds file
        out_dir (str): output directory
        cpu_number (int): number of CPU to use for multiprocessing

    Returns:
        set: union of reachable metabolites for all metabolic networks
    """
    # Run individual scopes of metabolic networks if any
    starttime = time.time()

    if len([
            name for name in os.listdir(sbmldir) if os.path.isfile(os.path.join(sbmldir, name))
            and utils.get_extension(os.path.join(sbmldir, name)).lower() in ["xml", "sbml"]
    ]) > 1:
        scope_dict, seeds_dict, scope_json, seeds_status_path = indiv_scope_run(sbmldir, seeds, out_dir, cpu_number)
        logger.info(f'\nIndividual scopes for all metabolic networks available in {scope_json}. The scopes have been filtered a way that if a seed is in a scope, it means the corresponding species is predicted to be able to produce it.')
        logger.info(f'\nInformation regarding the producibility of seeds, and the possible absence of seeds in some metabolic networks is stored in {seeds_status_path}.\n')        
        # Analyze the individual scopes results (json file)
        reachable_metabolites_union = analyze_indiv_scope(scope_dict, seeds_dict, seeds)
        # Compute the reverse iscopes (who produces each metabolite)
        reverse_scope_json, reverse_scope_tsv = reverse_scope(scope_dict, out_dir)
        logger.info(f"\nAnalysis of functional redundancy (producers of all metabolites) is computed as a dictionary in {reverse_scope_json} and as a matrix in {reverse_scope_tsv}.")
        logger.info("--- Indiv scopes runtime %.2f seconds ---\n" %
                    (time.time() - starttime))
        return reachable_metabolites_union
    else:
        logger.critical("The number of SBML files is <= 1 in the input directory")
        sys.exit(1)


def indiv_scope_run(sbml_dir, seeds, output_dir, cpu_number=1):
    """Run Menetools and analyse individual metabolic capabilities.
    
    Args:
        sbml_dir (str): directory of SBML files
        seeds (str): SBML seeds file
        output_dir (str): directory for results
        cpu_number (int): number of CPU to use for multiprocessing
    
    Returns:
        str: path to output file for scope from Menetools analysis
    """
    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#       Individual metabolic potentials       #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    menetools_dir = os.path.join(output_dir, 'indiv_scopes')
    indiv_scopes_path = os.path.join(menetools_dir, 'indiv_scopes.json')
    produced_seeds_path = os.path.join(menetools_dir, 'seeds_in_indiv_scopes.json')

    if not utils.is_valid_dir(menetools_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    all_files = [
        f for f in os.listdir(sbml_dir)
        if os.path.isfile(os.path.join(sbml_dir, f)) and utils.get_extension(
            os.path.join(sbml_dir, f)).lower() in ['xml', 'sbml']
    ]
    all_scopes = {}
    all_produced_seeds = {}
    all_absent_seeds = {}
    all_non_produced_seeds = {}
    multiprocessing_indiv_scopes = []
    for f in all_files:
        bname = utils.get_basename(f)
        sbml_path = os.path.join(sbml_dir, f)
        multiprocessing_indiv_scopes.append((sbml_path, bname, seeds))

    menescope_pool = Pool(cpu_number)
    results = menescope_pool.starmap(indiv_scope_on_species, multiprocessing_indiv_scopes)
    for result in results:
        error = result[0]
        if error is True:
            logger.critical('------------An error occurred during M2M run of Menetools, M2M will stop-------------')
            menescope_pool.close()
            menescope_pool.join()
            sys.exit(1)
        bname = result[1]
        menescope_results = result[2]
        all_scopes[bname] = menescope_results['scope']
        all_produced_seeds[bname] = menescope_results['produced_seeds']
        all_absent_seeds[bname] = menescope_results['absent_seeds']
        all_non_produced_seeds[bname] = menescope_results['non_produced_seeds']

    menescope_pool.close()
    menescope_pool.join()

    seeds_status = {}
    seeds_status['individually_producible_seeds'] = all_produced_seeds
    seeds_status['seeds_absent_in_metabolic_network'] = all_absent_seeds
    seeds_status['individually_non_producible_seeds'] = all_non_produced_seeds

    # Some seeds might be producible by some metabolic networks. By default, menescope include all seeds in the scope, regardless of their real producibility by the network. seed_status_dict holds this information.
    # We'll remove them here
    # remove from each scope the seeds that are not producible by the corresponding species. 
    for species in all_scopes:
        # non producible seeds are in seeds_status_dict['individually_non_producible_seeds'][species]
        all_scopes[species] = list(set(all_scopes[species]) - set(seeds_status['individually_non_producible_seeds'][species]))


    with open(indiv_scopes_path, 'w') as dumpfile:
        json.dump(all_scopes, dumpfile, indent=4)

    with open(produced_seeds_path, 'w') as dumpfile:
        json.dump(seeds_status, dumpfile, indent=4, sort_keys=True)

    return all_scopes, seeds_status, indiv_scopes_path, produced_seeds_path


def indiv_scope_on_species(sbml_path, bname, seeds_path):
    """Run Menetools and analyse individual metabolic capabilities on a sbml.

    Args:
        sbml_path (str): path to SBML file
        bname (str): name linked to SBML file
        seeds_path (str): path to SBML seeds file

    Returns:
        list: [boolean error, bname, dictionary containing menescope results]
    """
    error = False
    try:
        menescope_results = run_menescope(
            draft_sbml=sbml_path, seeds_sbml=seeds_path)
    except:
        traceback_str = traceback.format_exc()
        # Don't print the traceback if the error is linked to SystemExit as the error has been hanled by menetools.
        if 'SystemExit: 1' not in traceback_str:
            logger.critical(traceback_str)
        logger.critical('---------------Something went wrong running Menetools on ' + bname + '---------------')
        error = True
        menescope_results = None

    return [error, bname, menescope_results]


def analyze_indiv_scope(scope_dict, seeds_status_dict, seeds):
    """Analyze the output of Menescope, stored in two dictionaries
    
    Args:
        scope_dict (dict): output of all menescope runs
        seeds_status_dict (dict): production status of seeds in all menescope runs
        seeds (str): SBML seeds file

    Returns:
        set: union of all the individual scopes
    """
    scope_dict_set = {}
    individually_producible_seeds_set = {}

    for elem in scope_dict:
        scope_dict_set[elem] = set(scope_dict[elem])

    for elem in seeds_status_dict['individually_producible_seeds']:
        individually_producible_seeds_set[elem] = set(seeds_status_dict['individually_producible_seeds'][elem])

    try:
        seed_metabolites = readSBMLspecies_clyngor(seeds, 'seeds')
    except FileNotFoundError:
        logger.critical('File not found: '+ seeds)
        sys.exit(1)
    except etree.ParseError:
        logger.critical(f'Invalid syntax in SBML file: {seeds}')
        sys.exit(1)
    except:
        traceback_str = traceback.format_exc()
        # Don't print the traceback if the error is linked to SystemExit as the error has been handled by menetools.
        if 'SystemExit: 1' not in traceback_str:
            logger.critical(traceback_str)
        logger.critical('---------------Something went wrong running Menetools on " + seeds + "---------------')
        sys.exit(1)

    logger.info('%i metabolic models considered.' %(len(scope_dict_set)))
    intersection_scope = set.intersection(*list(scope_dict_set.values()))
    logger.info('\n' + str(len(intersection_scope)) + ' metabolites in core reachable by all organisms (intersection) \n')
    logger.info("\n".join(intersection_scope))

    union_scope = set.union(*list(scope_dict_set.values()))
    union_producible_seeds = set.union(*list(individually_producible_seeds_set.values()))
    logger.info('\n' + str(len(union_scope)) + ' metabolites reachable by individual organisms altogether (union), among which ' + str(len(union_producible_seeds)) + ' metabolites that are also part of the seeds (growth medium) \n')
    logger.info("\n".join(union_scope))
    len_scope = [len(scope_dict[elem]) for elem in scope_dict]
    logger.info('\nSummary:')
    logger.info('- intersection of scopes ' + str(len(intersection_scope)))
    logger.info('- union of scopes ' + str(len(union_scope)))
    logger.info('- max metabolites in scopes ' + str(max(len_scope)))
    logger.info('- min metabolites in scopes ' + str(min(len_scope)))
    logger.info('- average number of metabolites in scopes %.2f (+/- %.2f)' %
                (statistics.mean(len_scope), statistics.stdev(len_scope)))
    return union_scope


def reverse_scope(scope_dict, output_dir):
    """Reverse a scope dictionary by focusing on metabolite producers.

    Args:
        scope_dict (dict): dict of scope
        output_dir (str): path to output directory
    
    Returns:
        (str, str): paths to the JSON and TSV outputs
    """
    rev_indiv_scopes_json_path = os.path.join(*[output_dir, 'indiv_scopes', 'rev_iscope.json'])
    rev_indiv_scopes_tsv_path = os.path.join(*[output_dir, 'indiv_scopes', 'rev_iscope.tsv'])

    new_dic = {}
    for k,v in scope_dict.items():
        for x in v:
            new_dic.setdefault(x,[]).append(k)

    with open(rev_indiv_scopes_json_path, 'w') as g:
        json.dump(new_dic, g, indent=True, sort_keys=True)

    all_compounds = [compound for compound in new_dic]
    all_species = [species for species in scope_dict]

    # For each species get the possibility of production of each compounds.
    with open(rev_indiv_scopes_tsv_path, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['', *all_compounds])
        for species in all_species:
             csvwriter.writerow([species, *[1 if species in new_dic[compound] else 0 for compound in all_compounds]])

    return rev_indiv_scopes_json_path, rev_indiv_scopes_tsv_path
