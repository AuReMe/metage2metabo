from metage2metabo.m2m_workflow import run_workflow, recon, iscope, cscope, addedvalue, mincom, instance_community
from metage2metabo import sbml_management
import logging
import pkg_resources
import argparse
import pyasp
from metage2metabo import utils
import sys
import os
from shutil import copyfile
import re


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
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + LICENSE)

    # parent parser
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
        parents=[parent_parser_g, parent_parser_o, parent_parser_c],
        description=
        "Run metabolic network reconstruction for each annotated genome of the input directory, using Pathway Tools"
    )
    #TODO chose the SBML level for the output 2 or 3, default = 2
    #TODO create a seed command that creates a SBML file starting from a text file with compounds identifiers
    indivscope_parser = subparsers.add_parser(
        "iscope",
        help="individual scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_c
        ],
        description=
        "Compute individual scopes (reachable metabolites from seeds) for each metabolic network of the input directory"
    )
    comscope_parser = subparsers.add_parser(
        "cscope",
        help="community scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c
        ],
        description="Compute the community scope of all metabolic networks"
    )
    added_value_parser = subparsers.add_parser(
        "addedvalue",
        help="added value of microbiota's metabolism over individual's",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c
        ],
        description=
        "Compute metabolites that are reachable by the community/microbiota and not by individual organisms"
    )
    mincom_parser = subparsers.add_parser(
        "mincom",
        help="minimal communtity selection",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c
        ],
        description=
        "Select minimal-size community to make reachable a set of metabolites")
    mincom_parser.add_argument(
        "-t",
        "--targets",
        help="targets for metabolic analysis",
        required=True
        )
    seeds_parser = subparsers.add_parser(
        "seeds",
        help="creation of seeds SBML file",
        parents=[parent_parser_o],
        description=
        "Create a SBML file starting for a simple text file with metabolic compounds identifiers"
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
            parent_parser_c
        ],
        description=
        "Run the whole workflow: metabolic network reconstruction, individual and community scope analysis and community selection"
    )

    args = parser.parse_args()

    # test writing in out_directory if a subcommand is given else print version and help
    if args.cmd:
        if not utils.is_valid_dir(args.out):
            logger.critical("Impossible to access/create output directory")
            sys.exit(1)
    else:
        logger.info("m2m " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    #if modelhost is given as an arg: check the SBML level and turn it into 2 if needed
    if args.cmd in ["workflow", "mincom", "cscope", "addedvalue"] and args.modelhost:
        new_arg_modelhost = args.modelhost
        #TODO new_arg_modelhost = SbML2 version of args.modelhost if args.modelhost.getLevel!=2
    else:
        new_arg_modelhost = None


    # deal with given subcommand
    if args.cmd == "workflow":
        main_workflow(args.genomes, args.out, args.cpu, args.clean, args.seeds,
                      new_arg_modelhost)
    elif args.cmd in ["iscope","cscope","addedvalue","mincom"]:
        network_dir = check_sbml(args.networksdir, args.out)
        if args.cmd == "iscope":
            main_iscope(network_dir, args.seeds, args.out)
        elif args.cmd == "cscope":
            main_cscope(network_dir, args.seeds, args.out, new_arg_modelhost)
        elif args.cmd == "addedvalue":
            main_added_value(network_dir, args.seeds, args.out,
                             new_arg_modelhost)
        elif args.cmd == "mincom":
            main_mincom(network_dir, args.seeds, args.out, args.targets, new_arg_modelhost)
    elif args.cmd == "recon":
        main_recon(args.genomes, args.out, args.cpu, args.clean)
    elif args.cmd == "seeds":
        main_seeds(args.metabolites, args.out)


def main_workflow(*allargs):
    run_workflow(*allargs)


def main_recon(*allargs):
    pgdbdir, sbmldir = recon(*allargs)
    logger.info("PGDB created in " + pgdbdir)
    logger.info("SBML files created in " + sbmldir)


def main_iscope(*allargs):
    return iscope(*allargs)


def main_cscope(*allargs):
    comscope = cscope(*allargs)[1]
    logger.info(str(len(comscope)) + " metabolites reachable by the whole community/microbiota:")
    logger.info(', '.join(comscope))
    return comscope


def main_added_value(sbmldir, seeds, outdir, host):
    iscope_metabolites = main_iscope(sbmldir, seeds, outdir)
    logger.info(", ".join(iscope_metabolites))
    cscope_metabolites = main_cscope(sbmldir, seeds, outdir, host)
    newtargets = addedvalue(iscope_metabolites, cscope_metabolites)
    sbml_management.create_species_sbml(newtargets, outdir + "/community_analysis/targets.sbml")
    logger.info("Target file created with the addedvalue targets in: " +
                outdir + "/community_analysis/targets.sbml")


def main_mincom(sbmldir, seedsfiles, outdir, targets, host):
    #create instance
    instance = instance_community(sbmldir, seedsfiles, outdir, targets, host)
    #run mincom
    mincom(instance, outdir)


def main_seeds(metabolites_file, outdir):
    outfile = outdir + "/seeds.sbml"
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


def check_sbml(folder, outdir):
    all_files = [
        f for f in os.listdir(folder)
        if os.path.isfile(os.path.join(folder, f)) and utils.get_extension(
            os.path.join(folder, f)).lower() in ["xml", "sbml"]
    ]
    sbml_levels = {}
    make_new_sbmls = False
    for f in all_files:
        sbml_levels[f] = sbml_management.get_sbml_level(os.path.join(folder, f))
        if sbml_levels[f] != 2:
            make_new_sbmls = True
    if make_new_sbmls:
        sbml_dir = outdir + "/new_sbml/"
        logger.warning(
            "At least one SBML has not a suitable level for the tools. They will be transformed and created in "
            + sbml_dir + ". The others will be copied in this directory")
        if not utils.is_valid_dir(sbml_dir):
            logger.critical("Impossible to write in output directory")
            sys.exit(1)
        for f in all_files:
            if sbml_levels[f] != 2:
                #create level 2 SBML in sbml_dir
                sbml_management.sbml_to_sbml(
                    os.path.join(folder, f), os.path.join(sbml_dir, f), 2)
            else:
                #copy the original SBML in sbml_dir
                copyfile(os.path.join(folder, f), os.path.join(sbml_dir, f))
    else:
        sbml_dir = folder
    return sbml_dir


if __name__ == "__main__":
    main()