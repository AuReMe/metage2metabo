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
From metabolic network reconstruction with annotated genomes to metabolic capabilities screening to identify organisms of interest in a large microbiota.
"""
REQUIRES = """
Requirements here Pathway Tools installed and in $PATH, and NCBI Blast
"""

root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
logging.basicConfig(
        format='%(message)s', level=logging.INFO)  #%(name)s:%(message)s
logger = logging.getLogger(__name__)

# Check pyasp binaries.
pyasp_bin_path = pyasp.__path__[0] + '/bin/'
bin_check = []
for asp_bin in ['clasp', 'gringo3', 'gringo4']:
    bin_check.append(utils.is_valid_file(pyasp_bin_path + asp_bin))
if not all(bin_check):
    logger.critical("Error with pyasp installation, retry to install it:")
    logger.critical("pip install pyasp==1.4.3 --no-cache-dir --force-reinstall (add --user if you are not in a virtualenv)")
    sys.exit(1)


def main():
    """Run programm
    """
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
    parent_parser_l = argparse.ArgumentParser(add_help=False)
    parent_parser_l.add_argument(
        "-l",
        "--level",
        help="Level for SBML creation, 2 or 3",
        required=False,
        choices=[2,3],
        default=2)

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
            parent_parser_l
        ],
        description=
        "Run metabolic network reconstruction for each annotated genome of the input directory, using Pathway Tools"
    )
    indivscope_parser = subparsers.add_parser(
        "iscope",
        help="individual scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_c,
            parent_parser_q
        ],
        description=
        "Compute individual scopes (reachable metabolites from seeds) for each metabolic network of the input directory"
    )
    comscope_parser = subparsers.add_parser(
        "cscope",
        help="community scope computation",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c, parent_parser_q
        ],
        description="Compute the community scope of all metabolic networks")
    added_value_parser = subparsers.add_parser(
        "addedvalue",
        help="added value of microbiota's metabolism over individual's",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c, parent_parser_q
        ],
        description=
        "Compute metabolites that are reachable by the community/microbiota and not by individual organisms"
    )
    mincom_parser = subparsers.add_parser(
        "mincom",
        help="minimal communtity selection",
        parents=[
            parent_parser_n, parent_parser_s, parent_parser_o, parent_parser_m,
            parent_parser_c, parent_parser_q
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
        parents=[parent_parser_o, parent_parser_q],
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
            parent_parser_c, parent_parser_q
        ],
        description=
        "Run the whole workflow: metabolic network reconstruction, individual and community scope analysis and community selection"
    )

    args = parser.parse_args()
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
        logger.info("m2m " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    #if modelhost is given as an arg: check the SBML level and turn it into 2 if needed
    if args.cmd in ["workflow", "mincom", "cscope", "addedvalue"] and args.modelhost:
        new_arg_modelhost = check_sbml(args.modelhost, args.out, folder=False)
    else:
        new_arg_modelhost = None


    # deal with given subcommand
    if args.cmd == "workflow":
        main_workflow(args.genomes, args.out, args.cpu, args.clean, args.seeds,
                      new_arg_modelhost)
    elif args.cmd in ["iscope","cscope","addedvalue","mincom"]:
        if not os.path.isdir(args.networksdir):
            logger.critical(args.networksdir + " is not a correct directory path")
            sys.exit(1)
        network_dir = check_sbml(args.networksdir, args.out)
        if not utils.is_valid_file(args.seeds):
            logger.critical(args.seeds + " is not a correct filepath")
            sys.exit(1)
        else:
            if args.cmd == "iscope":
                main_iscope(network_dir, args.seeds, args.out)
            elif args.cmd == "cscope":
                main_cscope(network_dir, args.seeds, args.out, new_arg_modelhost)
            elif args.cmd == "addedvalue":
                main_added_value(network_dir, args.seeds, args.out,
                                new_arg_modelhost)
            elif args.cmd == "mincom":
                if not utils.is_valid_file(args.targets):
                    logger.critical(args.targets + " is not a correct filepath")
                    sys.exit(1)
                else:
                    main_mincom(network_dir, args.seeds, args.out, args.targets, new_arg_modelhost)
    elif args.cmd == "recon":
        main_recon(args.genomes, args.out, args.level, args.cpu, args.clean)
    elif args.cmd == "seeds":
        if not utils.is_valid_file(args.metabolites):
            logger.critical(args.metabolites + " is not a correct filepath")
            sys.exit(1)
        else:
            main_seeds(args.metabolites, args.out)


def main_workflow(*allargs):
    """Run main workflow
    """
    run_workflow(*allargs)


def main_recon(*allargs):
    """Run recon command
    """
    pgdbdir, sbmldir = recon(*allargs)
    logger.info("PGDB created in " + pgdbdir)
    logger.info("SBML files created in " + sbmldir)


def main_iscope(*allargs):
    """Run iscope command
    """
    return iscope(*allargs)


def main_cscope(*allargs):
    """Run cscope command
    """
    comscope = cscope(*allargs)[1]
    logger.info(str(len(comscope)) + " metabolites reachable by the whole community/microbiota:")
    logger.info(', '.join(comscope))
    return comscope


def main_added_value(sbmldir, seeds, outdir, host):
    """Run addedvalue command
    
    Args:
        sbmldir (str): SBML file directory
        seeds (str): SBML file for seeds
        outdir (str): results directory
        host (str): SBML file for host
    """
    iscope_metabolites = main_iscope(sbmldir, seeds, outdir)
    logger.info(", ".join(iscope_metabolites))
    cscope_metabolites = main_cscope(sbmldir, seeds, outdir, host)
    newtargets = addedvalue(iscope_metabolites, cscope_metabolites)
    sbml_management.create_species_sbml(newtargets, outdir + "/community_analysis/targets.sbml")
    logger.info("Target file created with the addedvalue targets in: " +
                outdir + "/community_analysis/targets.sbml")


def main_mincom(sbmldir, seedsfiles, outdir, targets, host):
    """Run mincom command
    
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
    mincom(instance, outdir)


def main_seeds(metabolites_file, outdir):
    """Run seeds command
    
    Args:
        metabolites_file (str): text file with metabolites IDs, one per line
        outdir (str): Results directory
    """
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


def check_sbml(inpt, outdir, folder = True):
    """Check whether one or several SBML level 3 files are in directory. If yes, convert them into a new directory and copy the SBML files that are correct into this same directory
    
    Args:
        inpt (str): SBML files directory
        outdir (str): Results directory
        folder (bool): Defaults to True. Change function behavior is input is a file or a folder
    
    Returns:
        str: filepath of file or directory, same as input if all SBMLs are level2
    """
    if folder:
        all_files = [
            f for f in os.listdir(inpt)
            if os.path.isfile(os.path.join(inpt, f)) and utils.get_extension(
                os.path.join(inpt, f)).lower() in ["xml", "sbml"]
        ]
        sbml_levels = {}
        make_new_sbmls = False
        for f in all_files:
            sbml_levels[f] = sbml_management.get_sbml_level(
                os.path.join(inpt, f))
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
                        os.path.join(inpt, f), os.path.join(sbml_dir, f), 2)
                else:
                    #copy the original SBML in sbml_dir
                    copyfile(os.path.join(inpt, f), os.path.join(sbml_dir, f))
        else:
            sbml_dir = inpt
        return sbml_dir
    else:
        if not utils.is_valid_file(inpt):
            logger.critical(inpt + " is not a correct filepath")
            sys.exit(1)
        else:
            sbml_level = sbml_management.get_sbml_level(inpt)
            if sbml_level != 2:
                newsbml = outdir + utils.get_basename(inpt) + "_lvl2.sbml"
                logger.warning(inpt + " was not in a suitable level for analysis. A converted file is created in " + newsbml)
                sbml_management.sbml_to_sbml(inpt, newsbml, 2)
            else:
                newsbml = inpt
            return newsbml


if __name__ == "__main__":
    main()