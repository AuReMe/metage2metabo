from padmet_utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet_utils.connection.sbmlGenerator import padmet_to_sbml, check
from padmet_utils.connection.sbml_to_sbml import from_sbml_to_sbml
from padmet_utils.connection.sbml_to_padmet import from_sbml_to_padmet
from padmet.utils.sbmlPlugin import convert_from_coded_id
from multiprocessing import Pool
from metage2metabo import utils
from libsbml import SBMLReader, writeSBMLToFile, SBMLDocument
import libsbml
import os
import logging
import sys


logger = logging.getLogger(__name__)


def get_sbml_level(sbml_file):
    """Get SBML Level of a file
    
    Args:
        sbml_file (str): SBML file
    
    Returns:
        int: SBML Level
    """
    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    return document.getLevel()


def create_species_sbml(metabolites, outputfile):
    """Create a SBML files with a list of species containing metabolites of the input set
    
    Args:
        metabolites (set): set of metabolites
        outputfile (str): SBML file to be written
    """
    document = libsbml.SBMLDocument(2, 1)
    model = document.createModel("metabolites")
    for compound in metabolites:
        name, stype, comp = convert_from_coded_id(compound)
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(compound), 'set species id')
        check(s.setName(name), 'set species name')
        check(s.setCompartment(comp), 'set species compartment')
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
    padmet = from_pgdb_to_padmet(
        pgdb_folder=species_pgdb_dir,
        extract_gene=True,
        no_orphan=noorphan_bool)

    padmet_to_sbml(padmet, species_sbml_file, sbml_lvl=sbml_level, verbose=False)

    sbml_check = utils.is_valid_path(species_sbml_file)
    return sbml_check


def pgdb_to_sbml(pgdb_dir, output_dir, noorphan_bool, sbml_level, cpu):
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

    logger.info("######### Creating SBML files #########")
    sbml_dir = output_dir + "/sbml"
    if not utils.is_valid_dir(sbml_dir):
        logger.critical("Impossible to access/create output directory")
        sys.exit(1)

    pgdb_to_sbml_pool = Pool(processes=cpu)

    multiprocess_data = []
    for species in os.listdir(pgdb_dir):
        multiprocess_data.append(
            [pgdb_dir + '/' + species,
            sbml_dir + '/' + species + '.sbml',
            sbml_level, noorphan_bool])

    sbml_checks = pgdb_to_sbml_pool.map(run_pgdb_to_sbml, multiprocess_data)

    pgdb_to_sbml_pool.close()
    pgdb_to_sbml_pool.join()

    if all(sbml_checks):
        return sbml_dir
    else:
        logger.critical("Error during padmet/sbml creation.")
        sys.exit(1)


def sbml_to_sbml(
        sbml_file,
        sbml_output_file,
        level_wanted,
        cpu=1
):
    """Transform a sbml into another level sbml

    Args:
        sbml_file (string): pathname to species sbml file
        sbml_output_file (string): pathname to output sbml file
        level_wanted (int): SBML level of the output SBML
        version (string): version of tthe database

    Returns:
        sbml_check (bool): Check if sbml file exists
    """
    from_sbml_to_sbml(sbml_file, sbml_output_file, level_wanted, cpu)

    sbml_check = utils.is_valid_path(sbml_output_file)
    return sbml_check