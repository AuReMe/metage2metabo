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
    """Get target from sbml

    Args:
        sbml_file (str): SBML file

    Returns:
        list: target
    """
    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()
    targets = [target.id for target in model.getListOfSpecies()]
    return targets

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
    """Create a SBML files with a list of species containing metabolites of the input set
    
    Args:
        metabolites (set): set of metabolites
        outputfile (str): SBML file to be written
    """
    document = libsbml.SBMLDocument(2, 1)
    model = document.createModel("metabolites")
    for compound in metabolites:
        compound = compound.strip('"')
        name, stype, comp = convert_from_coded_id(compound)
        s = model.createSpecies()
        sbmlGenerator.check(s, 'create species')
        sbmlGenerator.check(s.setId(compound), 'set species id')

        # Add name and compartment if found by padmet
        if name is not None:
            sbmlGenerator.check(s.setName(name), 'set species name')
        elif name is None:
            logger.warning("No name for " + compound)

        if comp is not None:
            sbmlGenerator.check(s.setCompartment(comp), 'set species compartment')
        elif comp is None:
            logger.warning("No compartment for " + compound)

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
