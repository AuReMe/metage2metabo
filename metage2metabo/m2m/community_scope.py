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
import os
import sys
import tempfile
import time

from metage2metabo import utils

from miscoto import run_scopes, run_instance

from shutil import copyfile

logger = logging.getLogger(__name__)
logging.getLogger("miscoto").setLevel(logging.CRITICAL)


def cscope(sbmldir, seeds, out_dir, targets_file=None, host=None):
    """Run community scope.
    
    Args:
        sbmldir (str): SBML files directory
        seeds (str): SBML file for seeds
        out_dir (str): output directory
        targets_file (str): targets file
        host (str, optional): Defaults to None. Host metabolic network (SBML)
    
    Returns:
        tuple: instance file (str) and community scope (set)
    """
    starttime = time.time()
    # Create instance for community analysis
    instance_com = instance_community(sbmldir, seeds, out_dir, targets_file, host)
    # Run community scope
    logger.info("Running whole-community metabolic scopes")
    community_reachable_metabolites = comm_scope_run(instance_com, out_dir, host)
    logger.info("--- Community scope runtime %.2f seconds ---\n" %
                (time.time() - starttime))
    return instance_com, community_reachable_metabolites


def instance_community(sbml_dir, seeds, output_dir, targets_file=None, host_mn=None):
    """Create ASP instance for community analysis.
    
    Args:
        sbml_dir (str): directory of symbionts SBML files
        seeds (str): seeds SBML file
        output_dir (str): directory for results
        targets_file (str): targets file
        host_mn (str): metabolic network file for host

    Returns:
        str: instance filepath
    """
    logger.info(
            "######### Creating metabolic instance for the whole community #########"
        )
    miscoto_dir = os.path.join(output_dir, 'community_analysis')
    if not utils.is_valid_dir(miscoto_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    # create tempfile
    fd, outputfile = tempfile.mkstemp(suffix='.lp', prefix='miscoto_', dir=miscoto_dir)

    instance_filepath = run_instance(
        bacteria_dir=sbml_dir,
        seeds_file=seeds,
        host_file=host_mn,
        targets_file=targets_file,
        output=outputfile)

    logger.info("Created instance in " + instance_filepath)
    return instance_filepath


def comm_scope_run(instance, output_dir, host_mn=None):
    """Run Miscoto_scope and analyse community metabolic capabilities
    
    Args:
        instance (str): instance filepath
        output_dir (str): directory for results
        host_mn (str): metabolic network file for host
    
    Returns:
        set: microbiota scope
    """
    miscoto_dir = os.path.join(output_dir, 'community_analysis')
    com_scopes_path = os.path.join(miscoto_dir, 'comm_scopes.json')

    if not utils.is_valid_dir(miscoto_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)
    microbiota_scope = run_scopes(lp_instance_file=instance)

    # Remove keys "host_prodtargets", "host_scope", "comhost_scope" and "host_unprodtargets" if there is no host:
    if host_mn is None:
        del microbiota_scope['host_prodtargets']
        del microbiota_scope['host_unprodtargets']
        del microbiota_scope['host_scope']
        del microbiota_scope['comhost_scope']

    with open(com_scopes_path, 'w') as dumpfile:
        json.dump(microbiota_scope, dumpfile, indent=4)

    logger.info('Community scopes for all metabolic networks available in ' +
                com_scopes_path)

    return set(microbiota_scope['com_scope'])
