import argparse
import logging
import os
import pkg_resources
import re
import sys
import tarfile
import time

from shutil import copyfile, which

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


from metage2metabo import sbml_management, utils
from metage2metabo.m2m_analysis import run_analysis_workflow, enumeration_analysis, stat_analysis, graph_analysis, powergraph_analysis, check_oog_jar_file

VERSION = pkg_resources.get_distribution("metage2metabo").version
LICENSE = """Copyright (C) Dyliss & Pleiade
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
metage2metabo is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.\n
"""
MESSAGE = """
Detection of key species among communities.
"""
REQUIRES = """
Requires: Oog jar file (http://www.biotec.tu-dresden.de/research/schroeder/powergraphs/download-command-line-tool.html) for powergraph visualization.
"""

root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
logging.basicConfig(
        format='%(message)s', level=logging.INFO)  #%(name)s:%(message)s
logger = logging.getLogger(__name__)

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
        help="Folder containing JSON file of single JSON file containing miscoto enumeration results",
        required=True,
        type = str)
    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        "-g",
        "--graph",
        metavar="GML_DIR_OR_FILE",
        help="Folder containing Graph GML file or single GML file",
        required=True)
    parent_parser_jar = argparse.ArgumentParser(add_help=False)
    parent_parser_jar.add_argument(
		"--oog",
		help="OOG jar file for powergraph",
		required=True,
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
        "Run miscoto enumeration on sbml species with seeds and targets"
    )
    stat_parser = subparsers.add_parser(
        "stats",
        help="statistics on key species",
        parents=[
            parent_parser_j, parent_parser_o, parent_parser_taxon, parent_parser_q
        ],
        description=
        "Compute statistics on key species in the community"
    )
    graph_parser = subparsers.add_parser(
        "graph",
        help="graph creation with enumeration solution",
        parents=[
            parent_parser_j, parent_parser_o, parent_parser_t, parent_parser_taxon, parent_parser_q
        ],
        description="Create the solution graph using the JSON from miscoto enumeration")
    powergraph_parser = subparsers.add_parser(
        "powergraph",
        help="powergraph creation and visualization",
        parents=[
            parent_parser_g, parent_parser_o, parent_parser_jar, parent_parser_q
        ],
        description=
        "Compress the GMl graph of solution and create a powergraph (bbl) and a svg of the graph"
    )
    wkf_parser = subparsers.add_parser(
        "workflow",
        help="whole workflow",
        parents=[
            parent_parser_s, parent_parser_n, parent_parser_t, parent_parser_m, parent_parser_o, parent_parser_jar,
            parent_parser_taxon, parent_parser_q
        ],
        description=
        "Run the whole workflow: miscoto enumeration, statistics on key species, graph on solution and powergraph creation"
    )

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # set up the logger
    if args.quiet:
        root_logger.setLevel(logging.CRITICAL)
    else:
        root_logger.setLevel(logging.INFO)

    # test writing in out_directory if a subcommand is given else print version and help
    if args.cmd:
        if not utils.is_valid_dir(args.out):
            logger.critical("Impossible to access/create output directory")
            sys.exit(1)
    else:
        logger.info("m2m_analysis " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

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


    # deal with given subcommand
    if args.cmd == "workflow":
        main_analysis_workflow(network_dir, args.targets, args.seeds, args.out, args.taxon,
                                args.oog, new_arg_modelhost)
    elif args.cmd == "enum":
        main_enumeration(network_dir, args.targets, args.seeds, args.out, new_arg_modelhost)
    elif args.cmd == "stats":
        main_stat(args.json, args.out, args.taxon)
    elif args.cmd == "graph":
        main_graph(args.json, args.targets, args.out, args.taxon)
    elif args.cmd == "powergraph":
        main_powergraph(args.graph, args.out, args.oog)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))


def main_analysis_workflow(*allargs):
    """Run main workflow
    """
    run_analysis_workflow(*allargs)


def main_enumeration(*allargs):
    """Run enumeration command
    """
    enumeration_analysis(*allargs)


def main_stat(*allargs):
    """Run stat command
    """
    stat_analysis(*allargs)


def main_graph(*allargs):
    """Run graph command
    """
    graph_analysis(*allargs)


def main_powergraph(*allargs):
    """Run powergraph command
    """
    powergraph_analysis(*allargs)


if __name__ == "__main__":
    main()