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
import os
import json
import sys
import subprocess
import time

from shutil import which

try:
    with open('powergrasp.cfg', 'w') as config_file:
        config_file.write("[powergrasp options]\n")
        config_file.write("SHOW_STORY = no\n")

    from powergrasp import compress_by_cc
    os.remove('powergrasp.cfg')
except ImportError:
    os.remove('powergrasp.cfg')
    raise ImportError('Requires powergrasp (https://github.com/Aluriak/PowerGrASP).')

try:
    import ete3
except ImportError:
    raise ImportError('Requires ete3 (https://github.com/etetoolkit/ete).')

from metage2metabo import utils
from metage2metabo import __version__ as VERSION

from metage2metabo.m2m_analysis.enumeration import enumeration_analysis
from metage2metabo.m2m_analysis.graph_compression import powergraph_analysis, check_oog_jar_file
from metage2metabo.m2m_analysis.solution_graph import graph_analysis
from metage2metabo.m2m_analysis.m2m_analysis_workflow import run_analysis_workflow

LICENSE = """Copyright (C) Dyliss & Pleiade
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
metage2metabo is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.\n
"""
MESSAGE = """
Detection of key species among communities.
"""
REQUIRES = """
Optional: Oog.jar file (https://github.com/AuReMe/metage2metabo/tree/main/external_dependencies) for powergraph svg creation.
"""

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Check ASP binaries.
if not which('clingo'):
    logger.critical('clingo is not in the Path, m2m_analysis can not work without it.')
    logger.critical('You can install with: pip install clyngor-with-clingo') 
    sys.exit(1)

def main():
    """Run programm
    """
    start_time = time.time()
    parser = argparse.ArgumentParser(
        "m2m_analysis",
        description=MESSAGE + " For specific help on each subcommand use: m2m_analysis {cmd} --help",
        epilog=REQUIRES, formatter_class=argparse.RawTextHelpFormatter
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
        help="cpu number for multi-process",
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
        help="Metabolic networks directory",
        required=True)
    parent_parser_t = argparse.ArgumentParser(add_help=False)
    parent_parser_t.add_argument(
		"-t",
		"--targets",
        metavar="TARGETS_DIR_OR_FILE",
		help="Folder containg sbml targets or single sbml file for metabolic analysis",
		required=True)
    parent_parser_m = argparse.ArgumentParser(add_help=False)
    parent_parser_m.add_argument(
        "-m",
        "--modelhost",
        help="Host metabolic model for community analysis",
        required=False,
        default=None)
    parent_parser_taxon = argparse.ArgumentParser(add_help=False)
    parent_parser_taxon.add_argument(
        "--taxon",
        help="Mpwt taxon file",
        required=False,
        default=None)
    parent_parser_j = argparse.ArgumentParser(add_help=False)
    parent_parser_j.add_argument(
        "-j",
        "--json",
        metavar="JSON_DIR_OR_FILE",
        help="Folder containing JSON files of single JSON file containing miscoto enumeration results",
        required=True,
        type = str)
    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        "-g",
        "--gml",
        metavar="GML_DIR_OR_FILE",
        help="Folder containing GML files of single GML file containing m2m_analysis graph results",
        required=True)
    parent_parser_jar = argparse.ArgumentParser(add_help=False)
    parent_parser_jar.add_argument(
		"--oog",
		help="OOG jar file for powergraph svg creation using Power Graph Command Line Tool",
		required=False,
		type=str)
    parent_parser_level = argparse.ArgumentParser(add_help=False)
    parent_parser_level.add_argument(
		"--level",
		help="Taxonomy level, must be: phylum, class, order, family, genus or species. By default, it is phylum.",
		required=False,
		type=str)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")
    enum_parser = subparsers.add_parser(
        "enum",
        help="enumeration using miscoto",
        parents=[
            parent_parser_s, parent_parser_n, parent_parser_t, parent_parser_m, parent_parser_o, parent_parser_q
        ],
        description=
        "Run miscoto enumeration on sbml species with seeds and targets",
        allow_abbrev=False
    )
    graph_parser = subparsers.add_parser(
        "graph",
        help="graph creation with enumeration solution",
        parents=[
            parent_parser_j, parent_parser_o, parent_parser_t, parent_parser_taxon, parent_parser_q,
            parent_parser_level
        ],
        description="Create the solution graph using the JSON from miscoto enumeration",
        allow_abbrev=False)
    powergraph_parser = subparsers.add_parser(
        "powergraph",
        help="powergraph creation and visualization",
        parents=[
            parent_parser_j, parent_parser_g, parent_parser_jar, parent_parser_q, parent_parser_taxon,
            parent_parser_level, parent_parser_o
        ],
        description=
        "Compress the GMl graph of solution and create a powergraph (bbl), a website format of the powergraph and a svg of the graph (if you use the --oog option)",
        allow_abbrev=False
    )
    wkf_parser = subparsers.add_parser(
        "workflow",
        help="whole workflow",
        parents=[
            parent_parser_s, parent_parser_n, parent_parser_t, parent_parser_m, parent_parser_o, parent_parser_jar,
            parent_parser_taxon, parent_parser_q, parent_parser_level
        ],
        description=
        "Run the whole workflow: miscoto enumeration, graph on solution and powergraph creation",
        allow_abbrev=False
    )

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # set up the logger
    if args.quiet:
        logger.setLevel(logging.CRITICAL)
    else:
        logger.setLevel(logging.INFO)

    # test writing in out_directory if a subcommand is given else print version and help
    if args.cmd:
        if not utils.is_valid_dir(args.out):
            logger.critical("Impossible to access/create output directory")
            sys.exit(1)
    else:
        logger.info("m2m_analysis " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    if "seeds" in args and args.seeds is not None:
        if not utils.is_valid_file(args.seeds):
            logger.critical('Error: ' + args.seeds + " is not a correct filepath")
            sys.exit(1)
    if "targets" in args and args.targets is not None:
        if not utils.is_valid_file(args.targets) and not utils.is_valid_dir(args.targets):
            logger.critical('Error: ' + args.targets + " is not a correct filepath or folder")
            sys.exit(1)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.out, f'm2m_analysis_{args.cmd}.log')
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

    # Check Oog.jar file
    if args.cmd in ["workflow", "powergraph"]:
        if args.oog:
            check_oog_jar_file(args.oog)

    #if modelhost is given as an arg: check the SBML level and turn it into 2 if needed
    if args.cmd in ["workflow", "enum"]:
        if not os.path.isdir(args.networksdir):
            logger.critical(args.networksdir + " is not a correct directory path")
            sys.exit(1)
        network_dir = args.networksdir

        if not utils.is_valid_file(args.seeds):
            logger.critical(args.seeds + " is not a correct filepath")
            sys.exit(1)
        if not utils.is_valid_file(args.targets) and not utils.is_valid_dir(args.targets):
            logger.critical(args.targets + " is not a correct filepath")
            sys.exit(1)
        if args.modelhost:
            new_arg_modelhost = args.modelhost
        else:
            new_arg_modelhost = None

    if args.cmd in ["workflow", "graph", "powergraph"]:
        if args.level:
            if args.level not in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
                logger.critical("Error with --level arugment, it must be one among: phylum, class, order, family, genus or species")
                sys.exit(1)
        if args.level is None:
            args.level = 'phylum'

    # deal with given subcommand
    if args.cmd == "workflow":
        main_analysis_workflow(network_dir, args.targets, args.seeds, args.out, args.taxon,
                                args.oog, new_arg_modelhost, args.level)
    elif args.cmd == "enum":
        main_enumeration(network_dir, args.targets, args.seeds, args.out, new_arg_modelhost)
    elif args.cmd == "graph":
        main_graph(args.json, args.targets, args.out, args.taxon, args.level)
    elif args.cmd == "powergraph":
        main_powergraph(args.json, args.gml, args.out, args.oog, args.taxon, args.level)

    duration = time.time() - start_time
    dict_args= vars(args)
    metadata_json_file = os.path.join(args.out, 'm2m_analysis_metadata.json')
    create_metadata(dict_args, duration, metadata_json_file)
    logger.info("--- Total runtime %.2f seconds ---" % (duration))


def main_analysis_workflow(*allargs):
    """Run main workflow
    """
    run_analysis_workflow(*allargs)


def main_enumeration(*allargs):
    """Run enumeration command
    """
    enumeration_analysis(*allargs)


def main_graph(*allargs):
    """Run graph command
    """
    graph_analysis(*allargs)


def main_powergraph(*allargs):
    """Run powergraph command
    """
    powergraph_analysis(*allargs)


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

    # If powergaph, get bubbletools and powergrasp versions.
    if dict_args['cmd'] == 'powergraph' or dict_args['cmd'] == 'workflow':
        if sys.version_info >= (3, 9):
            import importlib.metadata
            bubbletools_version = importlib.metadata.version("bubbletools")
            metadata['tool_dependencies']['python_package']['bubbletools'] = bubbletools_version
            powergrasp_version = importlib.metadata.version("powergrasp")
            metadata['tool_dependencies']['python_package']['powergrasp'] = powergrasp_version
        else:
            import pkg_resources
            bubbletools_version = pkg_resources.get_distribution('bubbletools').version
            metadata['tool_dependencies']['python_package']['bubbletools'] = bubbletools_version
            powergrasp_version = pkg_resources.get_distribution('powergrasp').version
            metadata['tool_dependencies']['python_package']['powergrasp'] = powergrasp_version

    # If graph, get networkx version.
    if dict_args['cmd'] == 'graph' or dict_args['cmd'] == 'workflow':
        from networkx import __version__ as networkx_version
        metadata['tool_dependencies']['python_package']['networkx'] = networkx_version

    # If taxonomy, get ete3 version.
    if dict_args['cmd'] == 'taxonomy' or dict_args['cmd'] == 'workflow':
        from ete3 import __version__ as ete3_version
        metadata['tool_dependencies']['python_package']['ete3'] = ete3_version

    metadata['duration'] = duration

    with open(metadata_json_file, 'w') as dumpfile:
        json.dump(metadata, dumpfile, indent=4)


if __name__ == "__main__":
    main()