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

import libsbml
import logging
import os
import sys

from libsbml import SBMLReader, writeSBMLToFile, SBMLDocument
from multiprocessing import Pool
from metage2metabo import utils
from padmet.utils.connection import pgdb_to_padmet, sbmlGenerator
from padmet.utils.sbmlPlugin import convert_from_coded_id

logger = logging.getLogger(__name__)


def get_compounds(sbml_file):
    """Get compound from sbml

    Args:
        sbml_file (str): SBML file

    Returns:
        list: compound
    """
    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()
    if model is None:
        logger.critical('SBML file "' + sbml_file + '" not well formatted. Is this file a SBML? Does it contains <model></model> tags?')
        sys.exit(1)
    compounds = [compound.id for compound in model.getListOfSpecies()]
    return compounds


def compare_seeds_and_targets(seedfile, targetfile):
    """Returns the intersection of the seeds and the targets

    Args:
        seedfile (str): path to seeds SBML file
        targetfile (str): path to targets SBML file

    Returns:
        set: intersection of seeds and targets
    """
    seeds = set(get_compounds(seedfile))
    targets = set(get_compounds(targetfile))

    return seeds.intersection(targets)


def create_species_sbml(metabolites, outputfile):
    """Create a SBML files with a list of species containing metabolites of the input set.
    Check if there are forbidden SBML characters in the metabolite IDs/ If yes, exit.
    
    Args:
        metabolites (set): set of metabolites
        outputfile (str): SBML file to be written
    """
    document = libsbml.SBMLDocument(2, 1)
    model = document.createModel("metabolites")
    forbidden_charlist = ['-', '|', '/', '(', ')',
        "'", '=', '#', '*', '.', ':', '!', '+', '[',
        ']', ',', ' ']
    forbidden_character_in_metabolites = None
    issue_trying_to_add_species = None
    for compound in metabolites:
        compound = compound.strip('"')
        name, stype, comp = convert_from_coded_id(compound)
        s = model.createSpecies()
        sbmlGenerator.check(s, 'create species')
        forbidden_characters_detacted = [char for char in forbidden_charlist if char in compound]
        if len(forbidden_characters_detacted) > 0:
            logger.warning("Forbidden character ({0}) in {1}. SBML creation will failed.".format(' '.join(forbidden_characters_detacted), compound))
            forbidden_character_in_metabolites = True
        try:
            sbmlGenerator.check(s.setId(compound), 'set species id')
        except:
            issue_trying_to_add_species = True
            logger.warning("Issue when trying to add compound {0}.".format(compound))

        if comp is not None:
            sbmlGenerator.check(s.setCompartment(comp), 'set species compartment')
        elif comp is None:
            logger.warning("No compartment for " + compound)

    if issue_trying_to_add_species is True and forbidden_character_in_metabolites is True:
        logger.warning("Forbidden character in compound ID, SBML creation will failed.")
        logger.warning("Modify the metabolic networks SBMl file by renaming these metabolites and removing the forbidden character.")
        sys.exit(1)
    if issue_trying_to_add_species is True and forbidden_character_in_metabolites is None:
        logger.warning("Issue when trying to add metabolite into SBML file, potential issue with SBML format.")
        logger.warning("Modify the metabolic networks SBMl file by renaming these metabolites and removing the forbidden character.")
        sys.exit(1)

    libsbml.writeSBMLToFile(document, outputfile)


def run_pgdb_to_sbml(species_multiprocess_data):
    """Turn PGDBs into SBML2 using multi-processing.

    Args:
        species_multiprocess_data (list): pathname to species pgdb dir, pathname to species sbml file

    Returns:
        sbml_check (bool): Check if sbml file exists
    """
    species_pgdb_dir = species_multiprocess_data[0]
    species_sbml_file = species_multiprocess_data[1]
    sbml_level = species_multiprocess_data[2]
    noorphan_bool = species_multiprocess_data[3]
    padmet_file_dir = species_multiprocess_data[4]

    padmet = pgdb_to_padmet.from_pgdb_to_padmet(
        pgdb_folder=species_pgdb_dir,
        extract_gene=True,
        no_orphan=noorphan_bool)

    if padmet_file_dir:
        padmet.generateFile(padmet_file_dir)

    sbmlGenerator.padmet_to_sbml(padmet, species_sbml_file, sbml_lvl=sbml_level, verbose=False)

    sbml_check = utils.is_valid_path(species_sbml_file)
    return sbml_check


def pgdb_to_sbml(pgdb_dir, output_dir, noorphan_bool, padmet_bool, sbml_level, cpu):
    """Turn Pathway Tools PGDBs into SBML2 files using Padmet
    
    Args:
        pgdb_dir (str): PGDB directory
        output_dir (str): results directory
        noorphan_bool (bool): ignores orphan reactions if True
        padmet_bool (bool): creates padmet files if True
        sbml_level (int): SBML level
        cpu (int): number of CPU for multi-process
    
    Returns:
        sbml_dir (str): SBML directory if successful
    """

    logger.info('######### Creating SBML files #########')
    sbml_dir = os.path.join(output_dir, 'sbml')
    padmet_dir = os.path.join(output_dir, 'padmet')

    if padmet_bool:
        if not utils.is_valid_dir(padmet_dir):
            logger.critical('Impossible to access/create output directory')
            sys.exit(1)
    if not utils.is_valid_dir(sbml_dir):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    pgdb_to_sbml_pool = Pool(processes=cpu)

    multiprocess_data = []
    for species in os.listdir(pgdb_dir):
        pgdb_species_path = os.path.join(pgdb_dir, species)
        sbml_species_path = os.path.join(sbml_dir, species + '.sbml')
        padmet_species_path = os.path.join(padmet_dir, species + '.padmet')
        if padmet_bool:
            multiprocess_data.append(
                [pgdb_species_path,
                sbml_species_path,
                sbml_level, noorphan_bool,
                padmet_species_path])
        else:
            multiprocess_data.append(
                [pgdb_species_path,
                sbml_species_path,
                sbml_level, noorphan_bool,
                padmet_bool])

    sbml_checks = pgdb_to_sbml_pool.map(run_pgdb_to_sbml, multiprocess_data)

    pgdb_to_sbml_pool.close()
    pgdb_to_sbml_pool.join()

    if all(sbml_checks):
        return sbml_dir
    else:
        logger.critical('Error during padmet/sbml creation.')
        sys.exit(1)
