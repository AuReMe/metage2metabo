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

import argparse
import logging
import json
import os
import pkg_resources
import re
import subprocess
import sys
import time
import traceback

from shutil import which

from metage2metabo import __version__ as VERSION
from metage2metabo.m2m.reconstruction import recon
from metage2metabo.m2m.individual_scope import iscope
from metage2metabo.m2m.community_scope import cscope, instance_community
from metage2metabo.m2m.community_addedvalue import addedvalue
from metage2metabo.m2m.minimal_community import mincom
from metage2metabo.m2m.m2m_workflow import run_workflow, metacom_analysis
from metage2metabo.sbml_management import get_compounds

from metage2metabo import sbml_management, utils

LICENSE = """Copyright (C) Dyliss & Pleiade\n
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.\n

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.\n

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>\n
"""
MESSAGE = """
From metabolic network reconstruction with annotated genomes to metabolic capabilities screening to identify organisms of interest in a large microbiota.
"""
REQUIRES = """
Requires: Pathway Tools installed and in $PATH, and NCBI Blast
"""

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Check ASP binaries.
if not which('clingo'):
    logger.critical('clingo is not in the Path, m2m can not work without it.')
    logger.critical('You can install with: pip install clyngor-with-clingo') 
    sys.exit(1)

# Check if clingo is accessible.
try:
    subprocess.call(['clingo', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except:
    logger.critical('Error with clingo.')
    logger.critical(traceback.format_exc())
    sys.exit(1)


def main():
    """Run programm.
    """
    start_time = time.time()
    parser = argparse.ArgumentParser(
        "m2m",
        description=MESSAGE + " For specific help on each subcommand use: m2m {cmd} --help",
        epilog=REQUIRES
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + LICENSE)

    # parent parser
    parent_parser_q = argparse.ArgumentParser(add_help=False)
    parent_parser_q.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        help="quiet mode",
        required=False,
        action="store_true",
        default=None,
    )
    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        "-c",
        "--cpu",
        help="cpu number for multiprocessing",
        required=False,
        type=int,
        default=1)
    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        "-o",
        "--out",
        dest="out",
        required=True,
        help="output directory path",
        metavar="OUPUT_DIR")

    parent_parser_no = argparse.ArgumentParser(add_help=False)
    parent_parser_no.add_argument(
        "--noorphan",
        help="use this option to ignore reactions without gene or protein association",
        required=False,
        action="store_true",
        default=False,
    )
    parent_parser_xml = argparse.ArgumentParser(add_help=False)
    parent_parser_xml.add_argument(
        "--pwt-xml",
        dest="pwt_xml",
        help="use this option to use Pathway Tools xml (incompatible with -p)",
        required=False,
        action="store_true",
        default=False,
    )
    parent_parser_s = argparse.ArgumentParser(add_help=False)
    parent_parser_s.add_argument(
        "-s",
        "--seeds",
        help="seeds (growth medium) for metabolic analysis",
        required=True)
    parent_parser_n = argparse.ArgumentParser(add_help=False)
    parent_parser_n.add_argument(
        "-n",
        "--networksdir",
        metavar="NETWORKS_DIR",
        help="metabolic networks directory",
        required=True)
    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        "-g", "--genomes", help="annotated genomes directory", required=True)
    parent_parser_cl = argparse.ArgumentParser(add_help=False)
    parent_parser_cl.add_argument(
        "--clean",
        help="clean PGDBs if already present",
        required=False,
        action="store_true",
        default=None)
    parent_parser_m = argparse.ArgumentParser(add_help=False)
    parent_parser_m.add_argument(
        "-m",
        "--modelhost",
        help="host metabolic model for community analysis",
        required=False,
        default=None)
    parent_parser_l = argparse.ArgumentParser(add_help=False)
    parent_parser_l.add_argument(
        "-l",
        "--level",
        help="Level for SBML creation, 2 or 3",
        required=False,
        type = int,
        choices=[2,3],
        default=2)
    parent_parser_p = argparse.ArgumentParser(add_help=False)
    parent_parser_p.add_argument(
        "-p",
        "--padmet",
        help="create padmet files",
        required=False,
        action="store_true",
        default=None)
    parent_parser_t_required = argparse.ArgumentParser(add_help=False)
    parent_parser_t_required.add_argument(
        "-t",
        "--targets",
        help="targets for metabolic analysis",
        required=True
        )
    parent_parser_t_optional = argparse.ArgumentParser(add_help=False)
    parent_parser_t_optional.add_argument(
        "-t",
        "--targets",
        help="Optional targets for metabolic analysis, if not used metage2metabo will use the addedvalue of the community",
        required=False
        )
    parent_parser_t_com_scope = argparse.ArgumentParser(add_help=False)
    parent_parser_t_com_scope.add_argument(
        "--target-com-scope",
        help="Instead of the addedvalue, use the community scope as targets for mincom.",
        required=False,
        action="store_true",
        default=None)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")
    ptools_parser = subparsers.add_parser(
        "recon",
        help="metabolic network reconstruction",
        parents=[
            parent_parser_g, parent_parser_o, parent_parser_c, parent_parser_q,
            parent_parser_l, parent_parser_no, parent_parser_p, parent_parser_cl,
            parent_parser_xml
        ],
        description=
        "Run metabolic network reconstruction for each annotated genome of the input directory, using Pathway Tools",
        allow_abbrev=False
    )
    indivscope_parser = subparsers.add_parser(
        "iscope",
        help="individual scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_q, parent_parser_c
        ],
        description=
        "Compute individual scopes (reachable metabolites from seeds) for each metabolic network of the input directory",
        allow_abbrev=False
    )
    comscope_parser = subparsers.add_parser(
        "cscope",
        help="community scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_q, parent_parser_t_optional
        ],
        description="Compute the community scope of all metabolic networks", allow_abbrev=False)
    added_value_parser = subparsers.add_parser(
        "addedvalue",
        help="added value of microbiota's metabolism over individual's",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_q
        ],
        description=
        "Compute metabolites that are reachable by the community/microbiota and not by individual organisms",
        allow_abbrev=False
    )
    mincom_parser = subparsers.add_parser(
        "mincom",
        help="minimal communtity selection",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_q, parent_parser_t_required
        ],
        description=
        "Select minimal-size community to make reachable a set of metabolites",
        allow_abbrev=False)
    seeds_parser = subparsers.add_parser(
        "seeds",
        help="creation of seeds SBML file",
        parents=[parent_parser_o, parent_parser_q],
        description=
        "Create a SBML file starting for a simple text file with metabolic compounds identifiers",
        allow_abbrev=False
    )
    seeds_parser.add_argument(
        "--metabolites",
        help=
        'metabolites file: one per line, encoded (XXX as in <species id="XXXX" .../> of SBML files)',
        required=True)
    wkf_parser = subparsers.add_parser(
        "workflow",
        help="whole workflow",
        parents=[
            parent_parser_g, parent_parser_s, parent_parser_m, parent_parser_o,
            parent_parser_c, parent_parser_q, parent_parser_no, parent_parser_p,
            parent_parser_t_optional, parent_parser_cl, parent_parser_xml,
            parent_parser_t_com_scope
        ],
        description=
        "Run the whole workflow: metabolic network reconstruction, individual and community scope analysis and community selection",
        allow_abbrev=False
    )
    metacom_parser = subparsers.add_parser(
        "metacom",
        help="whole metabolism community analysis",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_m, parent_parser_o,
            parent_parser_t_optional, parent_parser_q, parent_parser_c,
            parent_parser_t_com_scope
        ],
        description=
        "Run the whole metabolism community analysis: individual and community scope analysis and community selection",
        allow_abbrev=False
    )
    test_parser = subparsers.add_parser(
        "test",
        help="test on sample data from rumen experiments",
        parents=[
            parent_parser_q, parent_parser_c, parent_parser_o
        ],
        description="Test the whole workflow on a data sample", allow_abbrev=False)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # test writing in out_directory if a subcommand is given else print version and help
    if args.cmd:
        if not utils.is_valid_dir(args.out):
            logger.critical("Impossible to access/create output directory")
            sys.exit(1)
    else:
        logger.info("m2m " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.out, f'm2m_{args.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    if args.quiet:
        console_handler.setLevel(logging.WARNING)
    logger.addHandler(console_handler)

    #if modelhost is given as an arg: check the SBML level and turn it into 2 if needed
    if args.cmd in ["workflow", "metacom", "mincom", "cscope", "addedvalue"] and args.modelhost:
        new_arg_modelhost = args.modelhost
        logger.warning(f"\n A metabolic model is given for an host. The metabolite producibility of the community will display metabolites that can be produced by the hsot or the microbiome, *not including* the metabolites that the host can produce by itself. If this is not what you want to do, you can consider placing the host in the same directory as the symbionts, which will lead to a complete community scope. \n")
    else:
        new_arg_modelhost = None

    if args.cmd in ["workflow", "recon"]:
        if args.pwt_xml and args.padmet:
            logger.critical("-p and --pwt-xml are incompatible arguments")
            sys.exit(1)

    if "seeds" in args and args.seeds is not None:
        if not utils.is_valid_file(args.seeds):
            logger.critical(args.seeds + " is not a correct filepath")
            sys.exit(1)

    # deal with given subcommand
    if args.cmd == "workflow":
        main_workflow(args.genomes, args.out, args.cpu, args.clean, args.seeds,
                      args.noorphan, args.padmet, new_arg_modelhost, args.targets, args.pwt_xml,
                      args.target_com_scope)
    elif args.cmd in ["iscope", "cscope", "addedvalue", "mincom", "metacom"]:
        if not os.path.isdir(args.networksdir):
            logger.critical(args.networksdir + " is not a correct directory path")
            sys.exit(1)
        network_dir = args.networksdir
        if "targets" in args and args.targets is not None:
            if not utils.is_valid_file(args.targets):
                logger.critical(args.targets + " is not a correct filepath")
                sys.exit(1)
            # test if some targets are seeds
            itsct_seeds_targets = sbml_management.compare_seeds_and_targets(args.seeds, args.targets)
            if itsct_seeds_targets != set():
                logger.warning(f"\nWARNING: compounds {*list(itsct_seeds_targets),} are both in seeds and targets. As such, they will be considered already reachable during community selection and will be ignored. However, their producibility can be assessed in individual and community scopes. If a host is provided, the community scope computation differs slightly and the producibility of the compound will have to be checked in the output files: indiv_scopes/seeds_in_indiv_scopes.json and the key 'individually producible' in the file producibility_targets.json. \n")
        if args.cmd == "iscope":
            main_iscope(network_dir, args.seeds, args.out, args.cpu)
        elif args.cmd == "cscope":
            main_cscope(network_dir, args.seeds, args.out, args.targets, new_arg_modelhost)
        elif args.cmd == "addedvalue":
            main_added_value(network_dir, args.seeds, args.out,
                            new_arg_modelhost)
        elif args.cmd == "mincom":
            main_mincom(network_dir, args.seeds, args.out, args.targets, new_arg_modelhost)
        elif args.cmd == "metacom":
            main_metacom(network_dir, args.out, args.seeds, new_arg_modelhost, args.targets, args.cpu, args.target_com_scope)
    elif args.cmd == "recon":
        main_recon(args.genomes, args.out, args.noorphan, args.padmet, args.level, args.cpu,
                   args.clean, args.pwt_xml)
    elif args.cmd == "seeds":
        if not utils.is_valid_file(args.metabolites):
            logger.critical(args.metabolites + " is not a correct filepath")
            sys.exit(1)
        else:
            main_seeds(args.metabolites, args.out)
    elif args.cmd == 'test':
        main_test(args.out, args.cpu)

    duration = time.time() - start_time
    dict_args = vars(args)
    metadata_json_file = os.path.join(args.out, 'm2m_metadata.json')
    create_metadata(dict_args, duration, metadata_json_file)
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')


def main_workflow(*allargs):
    """Run main workflow.
    """
    run_workflow(*allargs)


def main_metacom(*allargs):
    """Run main workflow.
    """
    metacom_analysis(*allargs)


def main_recon(*allargs):
    """Run recon command.
    """
    pgdbdir, sbmldir, padmet_folder= recon(*allargs)
    logger.info("PGDB created in " + pgdbdir)
    if allargs[3] and padmet_folder:
        logger.info("PADMET files created in " + padmet_folder)

    logger.info("SBML files created in " + sbmldir)


def main_iscope(*allargs):
    """Run iscope command.
    """
    return iscope(*allargs)


def main_cscope(*allargs):
    """Run cscope command.
    """
    instance_com, comscope = cscope(*allargs)
    logger.info("\n" + str(len(comscope)) + " metabolites reachable by the whole community/microbiota: \n")
    logger.info('\n'.join(comscope))
    #delete intermediate file
    # Due to unstable behaviour of os.unlink on Windows, do not delete the file.
    # Refer to: https://github.com/python/cpython/issues/109608
    if sys.platform != 'win32':
        os.unlink(instance_com)
    return comscope


def main_added_value(sbmldir, seeds, outdir, host):
    """Run addedvalue command.
    
    Args:
        sbmldir (str): SBML file directory
        seeds (str): SBML file for seeds
        outdir (str): results directory
        host (str): SBML file for host
    """
    iscope_metabolites = main_iscope(sbmldir, seeds, outdir)
    #logger.info(", ".join(iscope_metabolites))
    cscope_metabolites = main_cscope(sbmldir, seeds, outdir, host)
    newtargets = addedvalue(iscope_metabolites, cscope_metabolites, outdir)
    if len(newtargets) > 0:
        targets_file_path = os.path.join(*[outdir, 'community_analysis', 'targets.sbml'])
        sbml_management.create_species_sbml(newtargets, targets_file_path)
        logger.info("Target file created with the addedvalue targets in: " +
                    targets_file_path)


def main_mincom(sbmldir, seedsfiles, outdir, targets, host):
    """Run mincom command.
    
    Args:
        sbmldir (str): SBML files directory
        seedsfiles (str): SBML file for seeds
        outdir (str): results directory
        targets (str): targets SBML file
        host (str): SBML file for host
    """
    #create instance
    instance = instance_community(sbmldir, seedsfiles, outdir, targets, host)
    #run mincom
    mincom(instance, seedsfiles, set(get_compounds(targets)), outdir)
    #delete intermediate file
    # Due to unstable behaviour of os.unlink on Windows, do not delete the file.
    # Refer to: https://github.com/python/cpython/issues/109608
    if sys.platform != 'win32':
        os.unlink(instance)


def main_seeds(metabolites_file, outdir):
    """Run seeds command.
    
    Args:
        metabolites_file (str): text file with metabolites IDs, one per line
        outdir (str): Results directory
    """
    outfile = os.path.join(outdir, 'seeds.sbml')
    with open(metabolites_file, "r") as f:
        rawdata = f.readlines()
    metabolites_set = set()
    for elem in rawdata:
        metabolites_set.add(elem.strip("\n"))
    for metabolite in metabolites_set:
        assert re.match(
            "^[A-Za-z0-9_-]*$", metabolite
        ) and not metabolite[0].isdigit(
        ), "Seed file is not in the correct format. Error with %s. Example of a correct ID is M_OXYGEN__45__MOLECULE_c. Rules = only numbers, letters or underscore in IDs, not starting with a number. One ID per line." %(metabolite)
    sbml_management.create_species_sbml(metabolites_set, outfile)
    logger.info("Seeds SBML file created in " + outfile)


def main_test(outdir, cpu):
    """Run test command.

    Args:
        outdir (str): directory containing the test data and the test output
        cpu (int): number of cpu to use (recommended: 2)
    """
    # Retrieve package path and path to test data.
    package_path = os.path.dirname(os.path.realpath(__file__))
    workflow_data_path = os.path.join(package_path, 'workflow_data')

    genome_file = os.path.join(workflow_data_path, 'workflow_genomes.tar.gz')
    seeds_sbml_file = os.path.join(workflow_data_path, 'seeds_workflow.sbml')

    logger.info("Uncompressing test data to " + outdir)
    utils.safe_tar_extract_all(genome_file, outdir)

    logger.info("Launching workflow on test data")
    input_genome = os.path.join(outdir, 'workflow_genomes')
    inp_dir=input_genome
    out_dir=outdir
    nb_cpu=cpu
    clean=None
    seeds=seeds_sbml_file
    noorphan_bool=None
    padmet_bool=True
    host_mn=None
    targets_file=None
    use_pwt_xml=False
    main_workflow(inp_dir, out_dir, nb_cpu,
                clean, seeds, noorphan_bool,
                padmet_bool, host_mn, targets_file,
                use_pwt_xml)


def create_metadata(dict_args, duration, metadata_json_file):
    """ Create metadata from args and package versions.

    Args:
        dict_args (dict): dict of args given to argparse
        duration (int): time of the run
        metadata_json_file (str): pat hto metadata output file
    """
    from miscoto import __version__ as miscoto_version
    from menetools import __version__ as menetools_version
    from metage2metabo import __version__ as m2m_version

    # Retrieve args given to m2m.
    metadata = {}
    metadata['m2m_args'] = dict_args

    # Get package version.
    metadata['tool_dependencies'] = {}
    metadata['tool_dependencies']['python_package'] = {}
    metadata['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata['tool_dependencies']['python_package']['Metage2Metabo'] = m2m_version
    metadata['tool_dependencies']['python_package']['MiSCoTo'] = miscoto_version
    metadata['tool_dependencies']['python_package']['MeneTools'] = menetools_version

    # Get clingo path and version.
    metadata['tool_dependencies']['clingo'] = which('clingo')
    response = subprocess.Popen(['clingo', '--version'], stdout=subprocess.PIPE, start_new_session=True, universal_newlines='')
    for line in response.stdout:
        str_line = str(line)
        if 'clingo version' in str_line:
            clingo_version = str_line.split('clingo version ')[1]
    metadata['tool_dependencies']['clingo_version'] = clingo_version

    # If recon, get mpwt, padmet and Pathway Tools versions.
    if dict_args['cmd'] == 'recon' or dict_args['cmd'] == 'workflow':
        from mpwt import __version__ as mpwt_version
        metadata['tool_dependencies']['python_package']['mpwt'] = mpwt_version
        if sys.version_info >= (3, 9):
            import importlib.metadata
            padmet_version = importlib.metadata.version("padmet")
        else:
            import pkg_resources
            padmet_version = pkg_resources.get_distribution('padmet').version

        metadata['tool_dependencies']['python_package']['padmet'] = padmet_version
        from mpwt.utils import get_ptools_version
        metadata['tool_dependencies']['Pathway-Tools'] = get_ptools_version()

    metadata['duration'] = duration

    with open(metadata_json_file, 'w') as dumpfile:
        json.dump(metadata, dumpfile, indent=4)


if __name__ == "__main__":
    main()
