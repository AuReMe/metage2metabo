#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

from shutil import copyfile

from metage2metabo import utils, sbml_management
from metage2metabo.m2m.reconstruction import recon
from metage2metabo.m2m.individual_scope import iscope
from metage2metabo.m2m.community_scope import cscope
from metage2metabo.m2m.community_addedvalue import addedvalue
from metage2metabo.m2m.minimal_community import mincom


logger = logging.getLogger(__name__)


def run_workflow(inp_dir, out_dir, nb_cpu, clean, seeds, noorphan_bool, padmet_bool, host_mn, targets_file, use_pwt_xml, target_com_scope=None):
    """Run the whole m2m workflow.
    
    Args:
        inp_dir (str): genomes directory
        out_dir (str): results directory
        nb_cpu (int): cpu number for multi-processing
        clean (bool): clean PGDB and re-run them
        seeds (str): seeds file
        noorphan_bool (bool): ignores orphan reactions if True
        padmet_bool (bool): creates padmet files if True
        host_mn (str): metabolic network file for host
        targets_file (str): targets file
        use_pwt_xml (bool): use Pathway Tools XML instead of creating them with padmet
        target_com_scope (bool): if True, will use all metabolties in com_scope as targets for minimal community predictions.
    """
    # METABOLIC NETWORK RECONSTRUCTION
    sbml_dir = recon(inp_dir, out_dir, noorphan_bool, padmet_bool, 2, nb_cpu, clean, use_pwt_xml)[1]

    # METABOLISM COMMUNITY ANALYSIS
    metacom_analysis(sbml_dir, out_dir, seeds, host_mn, targets_file, nb_cpu, target_com_scope)


def metacom_analysis(sbml_dir, out_dir, seeds, host_mn, targets_file, cpu_number=1, target_com_scope=None):
    """Run the metabolism community analysis part of m2m.

    Args:
        sbml_dir (str): sbml input directory
        out_dir (str): results directory
        seeds (str): seeds file
        host_mn (str): metabolic network file for host
        targets_file (str): targets file
        cpu_number (int): number of CPU to use for multiprocessing
        target_com_scope (bool): if True, will use all metabolties in com_scope as targets for minimal community predictions.
    """
    # INDIVIDUAL SCOPES
    union_targets_iscope = iscope(sbml_dir, seeds, out_dir, cpu_number)
    # COMMUNITY SCOPE
    instance_com, targets_cscope = cscope(sbml_dir, seeds, out_dir, targets_file, host_mn)
    # ADDED VALUE
    addedvalue_targets = addedvalue(union_targets_iscope, targets_cscope, out_dir)

    # If user gives a target file, check if the targets are in the addevalue
    if targets_file is not None:
        user_targets = set(sbml_management.get_compounds(targets_file))
        newtargets = user_targets
        individually_producible_targets = user_targets.intersection(union_targets_iscope)
        if len(individually_producible_targets) > 0:
            logger.info('\n The following ' + str(len(individually_producible_targets)) + " targets are individually reachable by at least one organism: \n")
            logger.info("\n".join(individually_producible_targets))
        commonly_producible_targets = user_targets.intersection(addedvalue_targets)
        if len(commonly_producible_targets) > 0:
            logger.info('\n The following ' + str(len(commonly_producible_targets)) + " targets are additionally reachable through putative cross-feeding events: \n")
            logger.info("\n".join(commonly_producible_targets))
        else:
            logger.info("Cross feeding interactions do not enable the producibility of additional targets")
    else:
        # If user has not specified to use the com_scope as targets, use the addedvalue.
        if target_com_scope is None:
            logger.info("\nUse the addedvalue as targets.")
            user_targets = None
            newtargets = addedvalue_targets
        else:
            logger.info("\nUse the community scope (all metabolites producible by the community) as targets.")
            user_targets = None
            newtargets = targets_cscope

    if len(newtargets) > 0:
        target_file_path = os.path.join(*[out_dir, 'community_analysis', 'targets.sbml'])
        if targets_file is not None:
            logger.info("\nTarget file created with the targets provided by the user in: " +
                        target_file_path)

        else:
            sbml_management.create_species_sbml(newtargets, target_file_path)
            if target_com_scope is None:
                logger_word = 'addedvalue'
            else:
                logger_word = 'community scope'
            logger.info("\nTarget file created with the {0} targets in: {1}".format(logger_word, target_file_path))

        sbml_management.create_species_sbml(newtargets, target_file_path)

        # Add these targets to the instance
        logger.info("Setting " + str(len(newtargets)) + " compounds as targets. \n")
        # if len(newtargets) != len(addedvalue_targets):
        #     logger.info("\n".join(newtargets))

        instance_w_targets = add_targets_to_instance(
            instance_com, out_dir,
            newtargets)
        # MINCOM
        mincom(instance_w_targets, seeds, newtargets, out_dir)
        # remove intermediate files
        # Due to unstable behaviour of os.unlink on Windows, do not delete the file.
        # Refer to: https://github.com/python/cpython/issues/109608
        if sys.platform != 'win32':
            os.unlink(instance_com)
        os.unlink(instance_w_targets)
    else:
        logger.info("No newly producible compounds, hence no community selection will be computed")
        os.unlink(instance_com)

    # Create targets producibility result file
    targets_producibility(out_dir, union_targets_iscope, targets_cscope, addedvalue_targets, user_targets, target_com_scope)


def add_targets_to_instance(instancefile, output_dir, target_set):
    """Add targets to the ASP community instance file.
    
    Args:
        instancefile (str): instance filepath
        output_dir (str): directory for results
        target_set (set): targets to be added
    
    Returns:
        str: new instance filepath
    """
    new_instance_file = os.path.join(*[output_dir, 'community_analysis', utils.get_basename(instancefile) + '__tgts.lp'])
    copyfile(instancefile, new_instance_file)

    with open(new_instance_file, 'a') as f:
        f.write('\n')
        for elem in target_set:
            f.write('target("' + elem + '").\n')

    return new_instance_file


def targets_producibility(m2m_out_dir, union_targets_iscope, targets_cscope, addedvalue_targets, user_targets=None, target_com_scope=None):
    """Create a json summarizing the producibility of the targets (either the addedvalue or the user provided targets)

    Args:
        m2m_out_dir (str): M2M results directory
        union_targets_iscope (list): targets producible by indiviual
        targets_cscope (list): targets producible by community
        addedvalue_targets (list): targets produbed by the community and not by individual
        user_targets (list): targets provided by the user
        target_com_scope (bool): if True, will use all metabolties in com_scope as targets for minimal community predictions.
    """
    prod_targets = {}

    if user_targets is not None:
        selected_targets = user_targets
        unproducible_targets = user_targets - union_targets_iscope - targets_cscope
        producible_targets = user_targets.intersection(union_targets_iscope.union(targets_cscope))
        indiv_producible = user_targets.intersection(union_targets_iscope)
    else:
        if target_com_scope is None:
            selected_targets = addedvalue_targets
            unproducible_targets = []
            producible_targets = addedvalue_targets
            indiv_producible = []
        else:
            selected_targets = targets_cscope
            unproducible_targets = []
            producible_targets = targets_cscope
            indiv_producible = []

    prod_targets['unproducible'] = list(unproducible_targets)
    prod_targets['producible'] = list(producible_targets)
    prod_targets['indiv_producible'] = list(indiv_producible)

    indiv_scopes_path = os.path.join(*[m2m_out_dir, 'indiv_scopes', 'indiv_scopes.json'])
    produced_seeds_path = os.path.join(*[m2m_out_dir, 'indiv_scopes', 'seeds_in_indiv_scopes.json'])
    comm_scopes_path = os.path.join(*[m2m_out_dir, 'community_analysis', 'comm_scopes.json'])
    reverse_cscope_path = os.path.join(*[m2m_out_dir, 'community_analysis', 'rev_cscope.json'])
    mincom_path = os.path.join(*[m2m_out_dir, 'community_analysis', 'mincom.json'])
    producibility_targets_path = os.path.join(m2m_out_dir, 'producibility_targets.json')

    if os.path.exists(indiv_scopes_path):
        prod_targets['individual_producers'] = {}
        with open(indiv_scopes_path) as json_data:
            producible_compounds = json.load(json_data)
        with open(produced_seeds_path) as json_data:
            producible_seeds = json.load(json_data)

        all_produced_seeds = []
        for species in producible_seeds:
            all_produced_seeds.extend(producible_seeds[species])
        all_produced_seeds = set(all_produced_seeds)

        for target in selected_targets:
            if target not in all_produced_seeds:
                species_producing_target = [species for species in producible_compounds if target in producible_compounds[species]]
                if species_producing_target != []:
                    prod_targets['individual_producers'][target] = species_producing_target
            else:
                species_producing_target = [species for species in producible_compounds if target in producible_seeds[species]]
                if species_producing_target != []:
                    prod_targets['individual_producers'][target] = species_producing_target

    if os.path.exists(comm_scopes_path):
        prod_targets['com_only_producers'] = {}
        if os.path.exists(reverse_cscope_path):
            with open(reverse_cscope_path) as json_data:
                rev_cscope = json.load(json_data)
            for target in selected_targets:
                if target in rev_cscope:
                    if target in prod_targets['individual_producers']:
                        only_com_producing_species = list(set(rev_cscope[target]) - set(prod_targets['individual_producers'][target]))
                    else:
                        only_com_producing_species = rev_cscope[target]
                    prod_targets['com_only_producers'][target] = only_com_producing_species
        else:
            with open(comm_scopes_path) as json_data:
                com_producible_compounds = json.load(json_data)
            if 'targets_producers' in com_producible_compounds and 'individual_producers' in com_producible_compounds:
                for target in selected_targets:
                    if target in com_producible_compounds['targets_producers']:
                        if target in prod_targets['individual_producers']:
                            only_com_producing_species = list(set(com_producible_compounds['targets_producers'][target]) - set(prod_targets['individual_producers'][target]))
                        else:
                            only_com_producing_species = com_producible_compounds['targets_producers'][target]
                        prod_targets['com_only_producers'][target] = only_com_producing_species

    if os.path.exists(mincom_path):
        with open(mincom_path) as json_data:
            mincom_producible_compounds = json.load(json_data)
        prod_targets['mincom_producible'] = mincom_producible_compounds['producible']
        prod_targets['key_species'] = mincom_producible_compounds['union_bacteria']
        prod_targets['alternative_symbionts'] = mincom_producible_compounds['alternative_symbionts']
        prod_targets['essential_symbionts'] = mincom_producible_compounds['essential_symbionts']
        prod_targets['mincom_optsol_producers'] = {}
        prod_targets['mincom_union_producers'] = {}
        prod_targets['mincom_inter_producers'] = {}
        for target in selected_targets:
            if target in mincom_producible_compounds['one_model_targetsproducers']:
                prod_targets['mincom_optsol_producers'][target] = mincom_producible_compounds['one_model_targetsproducers'][target]
            if target in mincom_producible_compounds['union_targetsproducers']:
                prod_targets['mincom_union_producers'][target] = mincom_producible_compounds['union_targetsproducers'][target]
            if target in mincom_producible_compounds['inter_targetsproducers']:
                prod_targets['mincom_inter_producers'][target] = mincom_producible_compounds['inter_targetsproducers'][target]

    with open(producibility_targets_path, 'w') as dumpfile:
        json.dump(prod_targets, dumpfile, indent=4)

    logger.info('Targets producibility are available at ' + producibility_targets_path)
