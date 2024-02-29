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

from metage2metabo import utils

logger = logging.getLogger(__name__)


def addedvalue(iscope_rm, cscope_rm, out_dir):
    """Compute the added value of considering interaction with microbiota metabolism rather than individual metabolisms.

    Args:
        iscope_rm (set): union of metabolites in all individual scopes
        cscope_rm (set): metabolites reachable by community/microbiota
        out_dir (str): output directory
    Returns:
        set: set of metabolites that can only be reached by a community
    """
    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#    Added-value of metabolic interactions    #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    # Community targets = what can be produced only if cooperation occurs between species
    newtargets = cscope_rm - iscope_rm
    logger.info("\nAdded value of cooperation over individual metabolism: " +
                str(len(newtargets)) + " newly reachable metabolites: \n")
    logger.info('\n'.join(newtargets))
    logger.info("\n")

    miscoto_dir = os.path.join(out_dir, 'community_analysis')
    addedvalue_json_path = os.path.join(miscoto_dir, 'addedvalue.json')

    if not utils.is_valid_dir(miscoto_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)
    dict_av = {'addedvalue': list(newtargets)}
    with open(addedvalue_json_path, 'w') as dumpfile:
        json.dump(dict_av, dumpfile, indent=4, default=lambda x: x.__dict__)
    logger.info(f'Added-value of cooperation written in {addedvalue_json_path}')

    return newtargets
