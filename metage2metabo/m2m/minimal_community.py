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

import json
import logging
import os
import sys
import time

from metage2metabo import utils
from metage2metabo.sbml_management import get_compounds

from miscoto import run_mincom

logger = logging.getLogger(__name__)
logging.getLogger("menetools").setLevel(logging.CRITICAL)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)
logging.getLogger("mpwt").setLevel(logging.INFO)


def mincom(instance_w_targets, seeds, targets, out_dir):
    """Compute minimal community selection and show analyses.
    
    Args:
        instance_w_targets (str): ASP instance filepath
        seeds (str): seeds filepath
        targets (str): targets set
        out_dir (str): results directory
    """
    starttime = time.time()

    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#         Minimal community selection         #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    miscoto_dir = os.path.join(out_dir, 'community_analysis')
    miscoto_mincom_path = os.path.join(miscoto_dir, 'mincom.json')

    # check if seeds are among the targets
    seeds = set(get_compounds(seeds))
    targets = set(targets)
    intersection = seeds.intersection(targets)
    if len(intersection) > 0:
        logger.warning(f'WARNING: The following seeds are among the targets: {intersection}. They will not be considered as targets during the computation of minimal communities: they will be considered as already reachable according to the network expansion definition.\n')    

    if not utils.is_valid_dir(miscoto_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)
    # Compute community selection
    logger.info('Running minimal community selection')
    all_results = compute_mincom(instance_w_targets, miscoto_dir)

    for key in all_results:
        all_results[key] = list(all_results[key])

    producible_targets = all_results['producible']
    unproducible_targets = all_results['still_unprod']
    logger.info('\nIn the initial and minimal communities ' + str(len(producible_targets)) + ' targets are producible and ' + str(len(unproducible_targets)) + ' remain unproducible.')
    logger.info('\n' + str(len(producible_targets)) + ' producible targets:') 
    logger.info('\n'.join(producible_targets))
    logger.info('\n' + str(len(unproducible_targets)) + ' still unproducible targets:') 
    logger.info('\n'.join(unproducible_targets))

    logger.info(f'\nMinimal communities are available in {miscoto_mincom_path} \n')
    # Give one solution
    one_sol_bact = []
    for bact in all_results['bacteria']:
        one_sol_bact.append(bact)
    logger.info('######### One minimal community #########')
    logger.info('# One minimal community enabling the producibility of the target metabolites given as inputs')
    logger.info('Minimal number of bacteria in communities => ' +
                str(len(one_sol_bact)) + '\n')
    logger.info("\n".join(one_sol_bact))
    # Give union of solutions
    union = all_results['union_bacteria']
    logger.info('######### Key species: Union of minimal communities #########')
    logger.info('# Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs')
    logger.info('Number of key species => ' +
                str(len(union)) + "\n")
    logger.info("\n".join(union))
    # Give intersection of solutions
    intersection = all_results['inter_bacteria']
    logger.info('######### Essential symbionts: Intersection of minimal communities #########')
    logger.info('# Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs')
    logger.info('Number of essential symbionts => ' +
                str(len(intersection)) + "\n")
    logger.info("\n".join(intersection))
    # Give key species, essential and alternative symbionts
    alternative_symbionts = list(set(union) - set(intersection))
    logger.info('######### Alternative symbionts: Difference between Union and Intersection #########')
    logger.info('# Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs')
    logger.info('Number of alternative symbionts => ' +
                str(len(alternative_symbionts)) + '\n')
    logger.info('\n'.join(alternative_symbionts))
    logger.info(
        '\n--- Mincom runtime %.2f seconds ---\n' % (time.time() - starttime))


def compute_mincom(instancefile, miscoto_dir):
    """Run minimal community selection and analysis.
    
    Args:
        instancefile (str): filepath to instance file
        miscoto_dir (str): directory with results

    Returns:
        dict: results of miscoto_mincom analysis
    """
    mincom_json_file = os.path.join(miscoto_dir, 'mincom.json')
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    results_dic = run_mincom(option="soup",
                            lp_instance_file=instancefile,
                            optsol=True,
                            union=True,
                            intersection=True,
                            output_json=mincom_json_file)
    return results_dic
