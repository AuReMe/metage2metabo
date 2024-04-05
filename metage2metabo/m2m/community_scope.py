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
import tempfile
import time
import csv

from metage2metabo import utils

from miscoto import run_focus, run_instance, run_scopes

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
    logger.info('\n###############################################')
    logger.info('#                                             #')
    logger.info('#    Metabolic potential of the community     #')
    logger.info('#                                             #')
    logger.info('###############################################\n')

    # Create instance for community analysis
    instance_com = instance_community(sbmldir, seeds, out_dir, targets_file, host)
    # Run community scope
    logger.info("Running whole-community metabolic scopes...")
    community_reachable_metabolites, contributions_of_microbes = comm_scope_run(instance_com, out_dir, host)
    # compute the reverse cscope
    contrib_microbes_path = os.path.join(*[out_dir, 'community_analysis', 'contributions_of_microbes.json'])
    # reverse the dict to have compounds as keys, and species as values
    reverse_contrib = {}
    if contributions_of_microbes is not None:
        for species in contributions_of_microbes:
            for compound in contributions_of_microbes[species]['produced_in_community']:
                if compound in reverse_contrib:
                    reverse_contrib[compound].append(species)
                else:
                    reverse_contrib[compound] = [species]
        # export the reverse cscope to json and tsv
        rev_cscopes_json_path, rev_cscopes_tsv_path = reverse_cscope(contributions_of_microbes, reverse_contrib, out_dir)
        logger.info('Reverse community scopes for all metabolic networks available in ' + rev_cscopes_json_path + ' and ' + rev_cscopes_tsv_path + '. They higlight the producibility of metabolites by species in the community.\n')
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

    logger.info("Created temporary instance file in " + instance_filepath)
    return instance_filepath


def comm_scope_run(instance, output_dir, host_mn=None):
    """Run Miscoto_scope and analyse community metabolic capabilities
    
    Args:
        instance (str): instance filepath
        output_dir (str): directory for results
        host_mn (str): metabolic network file for host
    
    Returns:
        set: microbiota scope
        dict: contribution of microbes to the scope
    """
    miscoto_dir = os.path.join(output_dir, 'community_analysis')
    com_scopes_path = os.path.join(miscoto_dir, 'comm_scopes.json')
    contrib_microbes_path = os.path.join(miscoto_dir, "contributions_of_microbes.json")

    if not utils.is_valid_dir(miscoto_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    # Remove keys "host_prodtargets", "host_scope", "comhost_scope" and "host_unprodtargets" if there is no host:
    if host_mn is not None:
        scopes_results = run_scopes(lp_instance_file=instance)
        com_scope_dict = {}
        com_scope_dict['com_scope'] = scopes_results['com_scope']
        com_scope_dict['host_prodtargets'] = scopes_results['host_prodtargets']
        com_scope_dict['host_scope'] = scopes_results['host_scope']
        com_scope_dict['comhost_scope'] = scopes_results['comhost_scope']
        com_scope_dict['host_unprodtargets'] = scopes_results['host_unprodtargets']

        contributions_of_microbes = None
        logger.info('The computation of the community scope with a host is of limited functionality. It will not highlight the contribution of each microbe to the community scope. Additionally, the producibility of seeds by the microbes will not be computed. Consider running the community scope without a host (i.e. all metabolic networks, including the host, in the same directory) to get the full functionality of the community scope. \n')
        with open(com_scopes_path, 'w') as dumpfile:
            json.dump(com_scope_dict, dumpfile, indent=4, sort_keys=True)
        logger.info(f'Community scope for all metabolic networks available in {com_scopes_path}.\n')
        microbiota_scope = set(com_scope_dict['com_scope'])
    else:
        contributions_of_microbes = run_focus(seeds_file = None, bacteria_dir = None, focus_bact=[], all_networks=True, lp_instance_file=instance)
        microbiota_scope = set()
        for bacteria in contributions_of_microbes:
            microbiota_scope.update(contributions_of_microbes[bacteria]['produced_in_community'])
        dict_comscope = {'com_scope': list(microbiota_scope)}
        with open(com_scopes_path, 'w') as dumpfile:
            json.dump(dict_comscope, dumpfile, indent=4, sort_keys=True)
        logger.info(f'Community scope for all metabolic networks available in {com_scopes_path}')
        with open(os.path.join(miscoto_dir, 'contributions_of_microbes.json'), 'w') as dumpfile:
            json.dump(contributions_of_microbes, dumpfile, indent=4, sort_keys=True)
        logger.info(f'Contributions of microbes to community scope available in {contrib_microbes_path}.\n')

        logger.info(f'\nNumber of metabolites producible in community: {len(microbiota_scope)}. \n')

    return microbiota_scope, contributions_of_microbes

def reverse_cscope(bact_contrib, reverse_dict, output_dir):
    """Reverse a scope dictionary by focusing on metabolite producers.

    Args:
        bact_contrib (dict): dict of bacteria contributions to community scope
        reverse_dict (dict): dict of metabolite producers in community
        output_dir (str): path to output directory
    
    Returns:
        (str, str): paths to the JSON and TSV outputs
    """
    rev_cscopes_json_path = os.path.join(*[output_dir, 'community_analysis', 'rev_cscope.json'])
    rev_cscopes_tsv_path = os.path.join(*[output_dir, 'community_analysis', 'rev_cscope.tsv'])

    with open(rev_cscopes_json_path, 'w') as g:
        json.dump(reverse_dict, g, indent=True, sort_keys=True)

    all_compounds = [compound for compound in reverse_dict]
    all_species = [species for species in bact_contrib]

    # For each species get the possibility of production of each compounds.
    with open(rev_cscopes_tsv_path, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['', *all_compounds])
        for species in all_species:
             csvwriter.writerow([species, *[1 if species in reverse_dict[compound] else 0 for compound in all_compounds]])

    return rev_cscopes_json_path, rev_cscopes_tsv_path
