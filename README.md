[![PyPI version](https://img.shields.io/pypi/v/metage2metabo.svg)](https://pypi.org/project/Metage2Metabo/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo.svg)](https://github.com/AuReMe/metage2metabo/blob/master/LICENSE) [![Actions Status](https://github.com/AuReMe/metage2metabo/workflows/Python%20package/badge.svg)](https://github.com/AuReMe/metage2metabo/actions) [![Documentation Status](https://readthedocs.org/projects/metage2metabo/badge/?version=latest)](https://metage2metabo.readthedocs.io/en/latest/?badge=latest) [![](https://img.shields.io/badge/doi-10.1101/803056-blueviolet.svg)](https://www.biorxiv.org/content/10.1101/803056v1)
# M2M - metage2metabo
Metage2metabo is a Python3 (Python >= 3.6) tool to perform graph-based metabolic analysis starting from annotated genomes (**reference genomes or metagenome-assembled genomes**). It uses *Pathway Tools* in a automatic and parallel way to **reconstruct metabolic networks** for a large number of genomes. The obtained metabolic networks are then **analyzed individually and collectively** in order to get the **added value of metabolic cooperation in microbiota over individual metabolism** and to **identify and screen interesting organisms** among all. 


m2m can be used as a whole workflow ( ``` m2m workflow ```, ``` m2m metacom ``` ) or steps can be performed individually ( ``` m2m recon ``` , ``` m2m iscope ``` , ``` m2m cscope ```, ``` m2m addedvalue ```, ``` m2m mincom ```, ``` m2m seeds ```).

## Table of contents
- [M2M - metage2metabo](#m2m---metage2metabo)
  - [Table of contents](#table-of-contents)
  - [General information about the modelling](#general-information-about-the-modelling)
  - [License](#license)
  - [Documentation](#documentation)
  - [Technologies](#technologies)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Installation with pip](#installation-with-pip)
    - [Availability on Docker and Singularity](#availability-on-docker-and-singularity)
  - [Features](#features)
    - [m2m recon](#m2m-recon)
    - [m2m iscope](#m2m-iscope)
    - [m2m cscope](#m2m-cscope)
    - [m2m addedvalue](#m2m-addedvalue)
    - [m2m mincom](#m2m-mincom)
    - [m2m workflow](#m2m-workflow)
    - [m2m metacom](#m2m-metacom)
    - [m2m seeds](#m2m-seeds)
    - [m2m test](#m2m-test)
  - [Analysis of the minimal solutions](#analysis-of-the-minimal-solutions)
  - [Release Notes](#release-notes)
  - [Additional features](#additional-features)
  - [Citation](#citation)
  - [Authors](#authors)
  - [Acknowledgement](#acknowledgement)

## General information about the modelling

M2M has two main dependencies for modelling metabolic networks: [MeneTools](https://github.com/cfrioux/MeneTools) and [Miscoto](https://github.com/cfrioux/miscoto). Accordingly metabolic models in M2M follow the producibility in metabolic networks as defined by the [network expansion](http://www.ncbi.nlm.nih.gov/pubmed/15712108) algorithm.
Mainly, two rules are followed:
* a *recursive rule*: the products of a reactions are producible if **all** reactants of this reaction are themselves producible
* an *initiation rule*: producibility is initiated by the presence of nutrients, called *seeds*. 

A metabolite that is producible from a set of nutrients is described as being "in the scope of the seeds".
The computation is made using logic solvers (Answer Set Programming). The present modelling ignores the stoichiometry of reactions (2A + B --> C is considered equivalent to A + B --> C), and is therefore suited to non-curated or draft metabolic networks, as the ones built using M2M with the PathoLogic software of [Pathway Tools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5036846/pdf/bbv079.pdf) handled by [Mpwt](https://github.com/AuReMe/mpwt). Many works have relied on network expansion to study organisms ([here](http://doi.wiley.com/10.1111/tpj.12627), [here](https://dx.plos.org/10.1371/journal.pcbi.1000049) or [there](http://dx.plos.org/10.1371/journal.pcbi.1005276)) and communities ([here](https://academic.oup.com/bioinformatics/article/34/17/i934/5093211), [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4786-7), or [here](https://www.ncbi.nlm.nih.gov/pubmed/18546499)). It has been [compared](http://www.ncbi.nlm.nih.gov/pubmed/19425125), [combined](https://www.cambridge.org/core/product/identifier/S1471068418000455/type/journal_article) to steady-state modelling (Flux Balance Analysis).

## License

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
* [padmet](https://github.com/AuReMe/padmet) to manage metabolic networks
* [menetools](https://github.com/cfrioux/MeneTools) to analyze individual metabolic capabilities using logic programming
* [miscoto](https://github.com/cfrioux/miscoto) to analyze collective metabolic capabilities and select communities within microbiota using logic programming

Also, m2m_analysis relies on other packages:
* [networkx](https://github.com/networkx/networkx) to create graph from miscoto results
* [ete3](https://github.com/etetoolkit/ete) to add taxonomy information on the graph if you used mpwt taxon file
* [powergrasp](https://github.com/Aluriak/PowerGrASP) to compress networkx graph

## Requirements

* [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) version 23.0 or higher (free for [academic users](https://biocyc.org/download-bundle.shtml)) is **required for m2m workflow and m2m recon**
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

* [Oog Power Graph Command line tool](http://www.biotec.tu-dresden.de/research/schroeder/powergraphs/download-command-line-tool.html) to create a svg file from the compressed graph at the end of m2m_analysis. This tool is a jar file, Java is needed to use it.

## Installation

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).
Tested with Python3.6

### Installation with pip

```
pip install Metage2Metabo
```

### Availability on Docker and Singularity

Due to Pathway-Tools license, Docker or Singularity images are not available publicly.

But you can create these images by using the Dockerfile and Singularity recipes available inside the recipes folder.
With these files, you can create container with Pathway-Tools and m2m.

More informations in the [Docker and Singularity Documentation](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-docker).

## Features

````
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
````

### m2m recon

````
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
````

### m2m iscope

````
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
````

### m2m cscope

````
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
````

### m2m addedvalue

````
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
````

### m2m mincom

````
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
````

### m2m workflow

````
usage: m2m workflow [-h] -g GENOMES -s SEEDS [-m MODELHOST] -o OUPUT_DIR
                    [-c CPU] [-q] [--noorphan] [-p] [-t TARGETS] [--clean]

Run the whole workflow: metabolic network reconstruction, individual and
community scope analysis and community selection

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMES, --genomes GENOMES
                        annotated genomes directory
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
  -p, --padmet          create padmet files
  -t TARGETS, --targets TARGETS
                        Optional targets for metabolic analysis, if not used
                        metage2metabo will use the addedvalue of the community
  --clean               clean PGDBs if already present
````

### m2m metacom

````
usage: m2m metacom [-h] -n NETWORKS_DIR -s SEEDS [-m MODELHOST] -o OUPUT_DIR
                   [-t TARGETS] [-q]

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
  -t TARGETS, --targets TARGETS
                        Optional targets for metabolic analysis, if not used
                        metage2metabo will use the addedvalue of the community
  -q, --quiet           quiet mode
````

### m2m seeds

````
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
````

### m2m test

````
usage: m2m test [-h] [-q] [-c CPU] -o OUPUT_DIR

Test the whole workflow on a data sample

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           quiet mode
  -c CPU, --cpu CPU     cpu number for multi-process
  -o OUPUT_DIR, --out OUPUT_DIR
                        output directory path
````

## Analysis of the minimal solutions

M2M performs a community minimization to find the union and intersection of the minimal communities. But it is possible to analyze all the minimal communities.
M2M has a second command-line, named m2m_analysis that performs this analysis. This method is slower than m2m as all sollutions are enumerated.
Then it creates a solutions graph and compresses it in a powergraph.

````
usage: m2m_analysis [-h] [-v] {enum,stats,graph,powergraph,workflow} ...

Detection of keystone species among communities.
 For specific help on each subcommand use: m2m_analysis {cmd} --help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid subcommands:

  {enum,stats,graph,powergraph,workflow}
    enum                enumeration using miscoto
    stats               statistics on keystone species
    graph               graph creation with enumeration solution
    powergraph          powergraph creation and visualization
    workflow            whole workflow

Requires: Oog jar file (http://www.biotec.tu-dresden.de/research/schroeder/powergraphs/download-command-line-tool.html) for powergraph visualization.
````

## Release Notes

Changes between version are listed on the [release page](https://github.com/AuReMe/metage2metabo/releases).

## Additional features

M2M relies on packages that can also be used independantly with more features:
* [mpwt](https://github.com/AuReMe/mpwt): command-line and multi-process solutions to run Pathway Tools. Suitable to multiple reconstruction, for example genomes of a microbiota
* [menetools](https://github.com/cfrioux/MeneTools): individual metabolic capabilities analysis using graph-based producibility criteria
* [miscoto](https://github.com/cfrioux/miscoto): community selection and metabolic screening in large-scal microbiotas, with or without taking a host into account

## Citation
Arnaud Belcour, Clémence Frioux, Meziane Aite, Anthony Bretaudeau, Anne Siegel (2019) Metage2Metabo: metabolic complementarity applied to genomes of large-scale microbiotas for the identification of keystone species. bioRxiv 803056; doi: [https://doi.org/10.1101/803056](https://doi.org/10.1101/803056)

## Authors
[Clémence Frioux](https://cfrioux.github.io/) and [Arnaud Belcour](https://arnaudbelcour.github.io/blog/), Univ Rennes, Inria, CNRS, IRISA, Rennes, France.

## Acknowledgement
People of Pathway Tools (SRI International) for their help integrating Pathway Tools with command line and multiprocessing in the [mpwt](https://github.com/AuReMe/mpwt) package, used in M2M.
