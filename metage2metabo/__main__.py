from metage2metabo.m2m_workflow import run_workflow
import logging
import pkg_resources
import argparse
import pyasp
from metage2metabo import utils
import sys

VERSION = pkg_resources.get_distribution("metage2metabo").version
LICENSE = """Copyright (C) Dyliss
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
m2m is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.\n
"""
MESSAGE = """
Description here TODO
"""
REQUIRES = """
Requirements here TODO
"""
logging.basicConfig(
        format='%(message)s', level=logging.INFO)  #%(name)s:%(message)s
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Check pyasp binaries.
pyasp_bin_path = pyasp.__path__[0] + '/bin/'
bin_check = []
for asp_bin in ['clasp', 'gringo3', 'gringo4']:
    bin_check.append(utils.is_valid_path(pyasp_bin_path + asp_bin))

if not all(bin_check):
    logger.critical("Error with pyasp installation, retry to install it:")
    logger.critical("pip install pyasp==1.4.3 --no-cache-dir --force-reinstall (add --user if you are not in a virtualenv)")
    sys.exit(1)


def main():

    parser = argparse.ArgumentParser(
        "m2m",
        description=MESSAGE + " For specific help on each subcommand use: m2m {cmd} --help",
        epilog=REQUIRES
    )
    parser.add_argument(
        "-c",
        "--cpu",
        help="cpu number for multi-process",
        required=False,
        type=int,
        default=1)
    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        default=False,
        help="do not print anything to standard output")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + LICENSE)

    # parent parser
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
        help="metabolic networks directory",
        required=True)
    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        "-g", "--genomes", help="annotated genomes directory", required=True)
    parent_parser_g.add_argument(
        "--clean",
        help="clean PGDBs if already present",
        required=False,
        action="store_true")
    parent_parser_m = argparse.ArgumentParser(add_help=False)
    parent_parser_m.add_argument(
        "-m",
        "--modelhost",
        help="host metabolic model for community analysis",
        required=False,
        default=None)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")
    ptools_parser = subparsers.add_parser(
        "recon",
        help="metabolic network reconstruction",
        parents=[parent_parser_g, parent_parser_o],
        description=
        "Run metabolic network reconstruction for each annotated genome of the input directory, using Pathway Tools"
    )
    indivscope_parser = subparsers.add_parser(
        "iscope",
        help="individual scope computation",
        parents=[parent_parser_n, parent_parser_s, parent_parser_o],
        description=
        "Compute individual scopes (reachable metabolites from seeds) for each metabolic network of the input directory"
    )
    comscope_parser = subparsers.add_parser(
        "cscope",
        help="community scope computation",
        parents=[parent_parser_n, parent_parser_s, parent_parser_o],
        description="Compute the community scope of all metabolic networks")

    added_value_parser = subparsers.add_parser(
        "addedvalue",
        help="added value of microbiota's metabolism over individual's",
        parents=[parent_parser_n, parent_parser_s, parent_parser_o],
        description=
        "Compute metabolites that are reachable by the community/microbiota and not by individual organisms"
    )
    mincom_parser = subparsers.add_parser(
        "mincom",
        help="minimal communtity selection",
        parents=[parent_parser_n, parent_parser_s, parent_parser_o,parent_parser_m],
        description=
        "Select minimal-size community to make reachable a set of metabolites")
    mincom_parser.add_argument(
        "-t",
        "--targets",
        help="targets for metabolic analysis",
        required=True
        )
    wkf_parser = subparsers.add_parser(
        "workflow",
        help="whole workflow",
        parents=[
            parent_parser_g, parent_parser_s, parent_parser_m, parent_parser_o
        ],
        description=
        "Run the whole workflow: metabolic network reconstruction, individual and community scope analysis and community selection"
    )


    args = parser.parse_args()

    ch = logging.StreamHandler()
    if args.quiet:
        ch.setLevel(logging.CRITICAL)
    else:
        ch.setLevel(logging.INFO)

    if args.cmd == "workflow":
        main_workflow(args.genomes, args.out, args.cpu, args.clean, args.seeds, args.modelhost)
    elif args.cmd == "recon":
        print("¯\_(ツ)_/¯ running recon")
    elif args.cmd == "iscope":
        print("¯\_(ツ)_/¯ running iscope")
    elif args.cmd == "cscope":
        print("¯\_(ツ)_/¯ running cscope")
    elif args.cmd == "addedvalue":
        print("¯\_(ツ)_/¯ running addedvalue")
    elif args.cmd == "mincom":
        print("¯\_(ツ)_/¯ running mincom")
    else:
        logger.info("m2m " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()


def main_workflow(*allargs):
    run_workflow(*allargs)


if __name__ == "__main__":
    main()