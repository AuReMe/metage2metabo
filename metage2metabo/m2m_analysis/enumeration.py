# Copyright (C) 2019-2021 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
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
import miscoto
import os
import sys
import time

from metage2metabo import utils

logger = logging.getLogger(__name__)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)


def enumeration(sbml_folder, target_file, seed_file, output_json, host_file):
    """Run miscoto enumeration on one target file

    Args:
        sbml_folder (str): sbml directory
        target_file (str): targets file
        seed_file (str): seeds file
        output_json (str): path to json output
        host_file (str): metabolic network file for host

    Returns:
        str: path to output json
    """
    results = miscoto.run_mincom(option="soup", bacteria_dir=sbml_folder,
        targets_file=target_file, seeds_file=seed_file,
        host_file=host_file, intersection=True,
        enumeration=True, union=True,
        optsol=True, output_json=output_json)

    # Give enumeration of solutions
    enumeration = str(len(results['enum_bacteria']))
    minimal_solution_size = str(len(results["bacteria"]))
    logger.info('######### Enumeration of minimal communities #########')
    logger.info(enumeration + ' minimal communities (each containing ' + minimal_solution_size + ' species) producing the target metabolites')
    # Give union of solutions
    union = results['union_bacteria']
    logger.info('######### Key species: Union of minimal communities #########')
    logger.info("# Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Key species = " +
                str(len(union)))
    logger.info("\n".join(union))
    # Give intersection of solutions
    intersection = results['inter_bacteria']
    logger.info('######### Essential symbionts: Intersection of minimal communities #########')
    logger.info("# Bacteria occurring in ALL minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Essential symbionts = " +
                str(len(intersection)))
    logger.info("\n".join(intersection))
    # Give key species, essential and alternative symbionts
    alternative_symbionts = list(set(union) - set(intersection))
    logger.info('######### Alternative symbionts: Difference between Union and Intersection #########')
    logger.info("# Bacteria occurring in at least one minimal community but not all minimal community enabling the producibility of the target metabolites given as inputs")
    logger.info("Alternative symbionts = " +
                str(len(alternative_symbionts)))
    logger.info("\n".join(alternative_symbionts))

    return output_json


def enumeration_analysis(sbml_folder, target_folder_file, seed_file, output_dir, host_file=None):
    """Run miscoto enumeration on input data

    Args:
        sbml_folder (str): sbml directory
        target_folder_file (str): targets file or folder containing multiple sbmls
        seed_file (str): seeds file
        output_dir (str): results directory
        host_file (str): metabolic network file for host

    Returns:
        dict: {target_filename_without_extension: json_output_path}
    """
    starttime = time.time()

    target_paths = utils.file_or_folder(target_folder_file)

    output_jsons = os.path.join(output_dir, 'json')
    if not utils.is_valid_dir(output_jsons):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    miscoto_jsons = {}
    for target_path in target_paths:
        logger.info('######### Enumeration of solution for: '+ target_path + ' #########')
        target_pathname = target_paths[target_path]
        output_json = os.path.join(output_jsons, target_path + '.json')
        if os.path.exists(output_json):
            logger.info('######### Enumeration has already been done for '+ target_path + ', it will not be launched again. #########')
        else:
            miscoto_json = enumeration(sbml_folder, target_pathname, seed_file, output_json, host_file)
            miscoto_jsons[target_path] = miscoto_json

    logger.info(
        "--- Enumeration runtime %.2f seconds ---\n" % (time.time() - starttime))

    return output_jsons
