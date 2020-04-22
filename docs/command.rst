========
Commands
========

Outputs and logs
-----------------

By default, ``m2m`` writes into the console (stdout) and into a log file located at the root of the results directory and named after the subcommand that wad executed. The option ``-q`` can be given to any ``m2m`` subcommand to write in the log file and write to stdout only the warnings, errors and critical issues, together with the path to the log file.

Features
========

    .. code::

        Copyright (C) Dyliss
        License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
        m2m is free software: you are free to change and redistribute it.
        There is NO WARRANTY, to the extent permitted by law.


        usage: m2m [-h] [-v]
                {recon,iscope,cscope,addedvalue,mincom,seeds,workflow,metacom,test}
                ...

        From metabolic network reconstruction with annotated genomes to metabolic
        capabilities screening to identify organisms of interest in a large
        microbiota. For specific help on each subcommand use: m2m {cmd} --help

        optional arguments:
        -h, --help            show this help message and exit
        -v, --version         show program's version number and exit

        subcommands:
        valid subcommands:

        {recon,iscope,cscope,addedvalue,mincom,seeds,workflow,metacom,test}
            recon               metabolic network reconstruction
            iscope              individual scope computation
            cscope              community scope computation
            addedvalue          added value of microbiota's metabolism over
                                individual's
            mincom              minimal communtity selection
            seeds               creation of seeds SBML file
            workflow            whole workflow
            metacom             whole metabolism community analysis
            test                test on sample data from rumen experiments

        Requires: Pathway Tools installed and in $PATH, and NCBI Blast


m2m recon
=========

    .. code::

        usage: m2m recon [-h] -g GENOMES -o OUPUT_DIR [-c CPU] [-q] [-l {2,3}]
                        [--noorphan] [-p] [--clean]

        Run metabolic network reconstruction for each annotated genome of the input
        directory, using Pathway Tools

        optional arguments:
            -h, --help            show this help message and exit
            -g GENOMES, --genomes GENOMES
                                    annotated genomes directory
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -c CPU, --cpu CPU     cpu number for multi-process
            -q, --quiet           quiet mode
            -l {2,3}, --level {2,3}
                                    Level for SBML creation, 2 or 3
            --noorphan            use this option to ignore reactions without gene or
                                    protein association
            -p, --padmet          create padmet files
            --clean               clean PGDBs if already present

m2m iscope
==========

    .. code::

        usage: m2m iscope [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-q]

        Compute individual scopes (reachable metabolites from seeds) for each
        metabolic network of the input directory

        optional arguments:
            -h, --help            show this help message and exit
            -n NETWORKS_DIR, --networksdir NETWORKS_DIR
                                    metabolic networks directory
            -s SEEDS, --seeds SEEDS
                                    seeds (growth medium) for metabolic analysis
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -q, --quiet           quiet mode

m2m cscope
==========

    .. code::

        usage: m2m cscope [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-m MODELHOST]
                    [-q]

        Compute the community scope of all metabolic networks

        optional arguments:
            -h, --help            show this help message and exit
            -n NETWORKS_DIR, --networksdir NETWORKS_DIR
                                    metabolic networks directory
            -s SEEDS, --seeds SEEDS
                                    seeds (growth medium) for metabolic analysis
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -m MODELHOST, --modelhost MODELHOST
                                    host metabolic model for community analysis
            -q, --quiet           quiet mode


m2m addedvalue
==============

    .. code::

        usage: m2m addedvalue [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR
                            [-m MODELHOST] [-q]

        Compute metabolites that are reachable by the community/microbiota and not by
        individual organisms

        optional arguments:
        -h, --help            show this help message and exit
        -n NETWORKS_DIR, --networksdir NETWORKS_DIR
                                metabolic networks directory
        -s SEEDS, --seeds SEEDS
                                seeds (growth medium) for metabolic analysis
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory path
        -m MODELHOST, --modelhost MODELHOST
                                host metabolic model for community analysis
        -q, --quiet           quiet mode


m2m mincom
==========

    .. code::

        usage: m2m mincom [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-m MODELHOST]
                    [-q] -t TARGETS

        Select minimal-size community to make reachable a set of metabolites

        optional arguments:
            -h, --help            show this help message and exit
            -n NETWORKS_DIR, --networksdir NETWORKS_DIR
                                    metabolic networks directory
            -s SEEDS, --seeds SEEDS
                                    seeds (growth medium) for metabolic analysis
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -m MODELHOST, --modelhost MODELHOST
                                    host metabolic model for community analysis
            -q, --quiet           quiet mode
            -t TARGETS, --targets TARGETS
                                targets for metabolic analysis

m2m workflow
============

    .. code::

        usage: m2m workflow [-h] -g GENOMES [--clean] -s SEEDS [-m MODELHOST] -o
                            OUPUT_DIR [-c CPU] [-q] [--noorphan]

        Run the whole workflow: metabolic network reconstruction, individual and
        community scope analysis and community selection

        optional arguments:
            -h, --help            show this help message and exit
            -g GENOMES, --genomes GENOMES
                                    annotated genomes directory
            --clean               clean PGDBs if already present
            -s SEEDS, --seeds SEEDS
                                    seeds (growth medium) for metabolic analysis
            -m MODELHOST, --modelhost MODELHOST
                                    host metabolic model for community analysis
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -c CPU, --cpu CPU     cpu number for multi-process
            -q, --quiet           quiet mode
            --noorphan            use this option to ignore reactions without gene or
                                    protein association


m2m metacom
===========

    .. code::

        usage: m2m metacom [-h] -n NETWORKS_DIR -s SEEDS [-m MODELHOST] -o OUPUT_DIR
                        [-q]

        Run the whole metabolism community analysis: individual and community scope
        analysis and community selection

        optional arguments:
        -h, --help            show this help message and exit
        -n NETWORKS_DIR, --networksdir NETWORKS_DIR
                                metabolic networks directory
        -s SEEDS, --seeds SEEDS
                                seeds (growth medium) for metabolic analysis
        -m MODELHOST, --modelhost MODELHOST
                                host metabolic model for community analysis
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory path
        -q, --quiet           quiet mode


m2m seeds
=========

    .. code::

        usage: m2m seeds [-h] -o OUPUT_DIR [-q] --metabolites METABOLITES

        Create a SBML file starting for a simple text file with metabolic compounds
        identifiers

        optional arguments:
            -h, --help            show this help message and exit
            -o OUPUT_DIR, --out OUPUT_DIR
                                    output directory path
            -q, --quiet           quiet mode
            --metabolites METABOLITES
                                    metabolites file: one per line, encoded (XXX as in
                                    <species id="XXXX" .../> of SBML files)
