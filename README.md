[![PyPI version](https://img.shields.io/pypi/v/metage2metabo.svg)](https://pypi.org/project/Metage2Metabo/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo.svg)](https://github.com/AuReMe/metage2metabo/blob/master/LICENSE) [![Build Status](https://travis-ci.org/AuReMe/metage2metabo.svg?branch=master)](https://travis-ci.org/AuReMe/metage2metabo)
# M2M - metage2metabo
Metage2metabo is a Python3 tool to perform graph-based metabolic analysis starting from annotated genomes (**reference genomes or metagenome-assembled genomes**). It uses *Pathway Tools* in a automatic and parallel way to **reconstruct metabolic networks** for a large number of genomes. The obtained metabolic networks are then **analyzed individually and collectively** in order to get the **added value of metabolic cooperation in microbiota over individual metabolism** and to **identify and screen interesting organisms** among all. 

m2m can be used as a whole workflow ( ``` m2m workflow ``` ) or steps can be performed individually ( ``` m2m recon ``` , ``` m2m iscope ``` , ``` m2m cscope ```, ``` m2m addedvalue ```, ``` m2m mincom ```, ``` m2m seeds ```).

## Table of contents
- [M2M - metage2metabo](#m2m---metage2metabo)
  - [Table of contents](#table-of-contents)
  - [General information](#general-information)
  - [Documentation](#documentation)
  - [Technologies](#technologies)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Installation with pip](#installation-with-pip)
    - [Availability on Docker and Singularity](#availability-on-docker-and-singularity)
  - [Features`](#features)
    - [m2m recon](#m2m-recon)
    - [m2m iscope](#m2m-iscope)
    - [m2m cscope](#m2m-cscope)
    - [m2m addedvalue](#m2m-addedvalue)
    - [m2m mincom](#m2m-mincom)
    - [m2m workflow](#m2m-workflow)
    - [m2m seeds](#m2m-seeds)
  - [Additional features](#additional-features)
  - [Authors](#authors)
  - [Acknowledgement](#acknowledgement)

## General information

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/AuReMe/metage2metabo/blob/master/LICENSE) file for details.

## Documentation

A more detailled documentation is available at: [https://metage2metabo.readthedocs.io](https://metage2metabo.readthedocs.io/en/latest/).

## Technologies

Python 3 (Python 3.6 is tested). M2M uses a certain number of Python dependencies. An example of all these dependencies working for Ubuntu 18.04 is available in [requirements.txt](https://github.com/AuReMe/metage2metabo/blob/master/requirements.txt).
They can be installed with:
````sh
pip install -r requirements.txt --no-cache-dir
````
In particular, m2m relies on:
* [mpwt](https://github.com/AuReMe/mpwt) to automatize metabolic network reconstruction with Pathway Tools
* [padmet](https://github.com/AuReMe/padmet) and [padmet-utils](https://github.com/AuReMe/padmet-utils) to manage metabolic networks
* [menetools](https://github.com/cfrioux/MeneTools) to analyze individual metabolic capabilities using logic programming
* [miscoto](https://github.com/cfrioux/miscoto) to analyze collective metabolic capabilities and select communities within microbiota using logic programming

## Requirements

* [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) version 22.5 or higher (free for [academic users](https://biocyc.org/download-bundle.shtml))
    * Pathway Tools requirements
        * **Linux**: Gnome terminal and Libxm4
        ```sh
        apt-get update && apt-get install gnome-terminal libxm4
         ```
        * **All OS**: [NCBI Blast](https://www.ncbi.nlm.nih.gov/books/NBK279671/) and a ncbirc file in user's home directory
            * Install with apt-get
            ```sh
            apt-get update && apt-get install gnome-terminal libxm4 ncbi-blast+ 
            echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc
            ```
            * Install with a [dmg installer](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) on MacOS
    * Pathway Tools install
        * **Linux**
        ```sh
        chmod +x ./pathway-tools-22.5-linux-64-tier1-install 
        ./pathway-tools-22.5-linux-64-tier1-install 
        ```
        and follow the instructions during the interactive install

        *For a silent install*: ```./pathway-tools-22.5-linux-64-tier1-install --InstallDir your/install/directory/pathway-tools --PTOOLS_LOCAL_PATH your/chosen/directory/for/data/ptools --InstallDesktopShortcuts 0 --mode unattended```
        * **MacOS**

        Dmg installer with a graphical interface.

        * **Warning** 
    
        /!\ For all OS, Pathway Tools must be in ```$PATH```. 
        On Linux and MacOS: ```export PATH=$PATH:your/install/directory/pathway-tools```. 
        Consider adding Pathway Tools in ```$PATH``` permanently by running
        ````sh
        echo 'export PATH="$PATH:your/install/directory/pathway-tools:"' >> ~/.bashrc
        ````

## Installation

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).
Tested with Python3.6

### Installation with pip

```pip install Metage2Metabo```

### Availability on Docker and Singularity

```details ```

## Features
````
Copyright (C) Dyliss
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
m2m is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.


usage: m2m [-h] [-v]
           {recon,iscope,cscope,addedvalue,mincom,seeds,workflow} ...

From metabolic network reconstruction with annotated genomes to metabolic
capabilities screening to identify organisms of interest in a large
microbiota. For specific help on each subcommand use: m2m {cmd} --help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid subcommands:

  {recon,iscope,cscope,addedvalue,mincom,seeds,workflow}
    recon               metabolic network reconstruction
    iscope              individual scope computation
    cscope              community scope computation
    addedvalue          added value of microbiota's metabolism over
                        individual's
    mincom              minimal communtity selection
    seeds               creation of seeds SBML file
    workflow            whole workflow

Requirements here Pathway Tools installed and in $PATH, and NCBI Blast
````

### m2m recon

````
usage: m2m recon [-h] -g GENOMES [--clean] -o OUPUT_DIR [-c CPU]

Run metabolic network reconstruction for each annotated genome of the input
directory, using Pathway Tools

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMES, --genomes GENOMES
                        annotated genomes directory
  --clean               clean PGDBs if already present
  -o OUPUT_DIR, --out OUPUT_DIR
                        output directory path
  -c CPU, --cpu CPU     cpu number for multi-process
  ````

### m2m iscope

````
usage: m2m iscope [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-c CPU]

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
  -c CPU, --cpu CPU     cpu number for multi-process
  ````

### m2m cscope

````
usage: m2m cscope [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-m MODELHOST]
                  [-c CPU]

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
  -c CPU, --cpu CPU     cpu number for multi-process
  ````

### m2m addedvalue

````
usage: m2m addedvalue [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR
                      [-m MODELHOST] [-c CPU]

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
  -c CPU, --cpu CPU     cpu number for multi-process
  ````

### m2m mincom

````
usage: m2m mincom [-h] -n NETWORKS_DIR -s SEEDS -o OUPUT_DIR [-m MODELHOST]
                  [-c CPU] -t TARGETS

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
  -c CPU, --cpu CPU     cpu number for multi-process
  -t TARGETS, --targets TARGETS
                        targets for metabolic analysis
````

### m2m workflow

````
usage: m2m workflow [-h] -g GENOMES [--clean] -s SEEDS [-m MODELHOST] -o
                    OUPUT_DIR [-c CPU]

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
  ````

### m2m seeds

````
usage: m2m seeds [-h] -o OUPUT_DIR --metabolites METABOLITES

Create a SBML file starting for a simple text file with metabolic compounds
identifiers

optional arguments:
  -h, --help            show this help message and exit
  -o OUPUT_DIR, --out OUPUT_DIR
                        output directory path
  --metabolites METABOLITES
                        metabolites file: one per line, encoded (XXX as in
                        <species id="XXXX" .../> of SBML files
````

## Additional features

M2M relies on packages that can also be used independantly with more features:
* [mpwt](https://github.com/AuReMe/mpwt): command-line and multi-process solutions to run Pathway Tools. Suitable to multiple reconstruction, for example genomes of a microbiota
* [menetools](https://github.com/cfrioux/MeneTools): individual metabolic capabilities analysis using graph-based producibility criteria
* [miscoto](https://github.com/cfrioux/miscoto): community selection and metabolic screening in large-scal microbiotas, with or without taking a host into account

## Authors
[Clémence Frioux](https://cfrioux.github.io/) and [Arnaud Belcour](https://arnaudbelcour.github.io/blog/), Univ Rennes, Inria, CNRS, IRISA, Rennes, France.

## Acknowledgement
People of Pathway Tools (SRI International) for their help integrating Pathway Tools with command line and multi process in the [mpwt](https://github.com/AuReMe/mpwt) package, used in M2M.
