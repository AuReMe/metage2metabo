[![PyPI version](https://img.shields.io/pypi/v/metage2metabo.svg)](https://pypi.org/project/Metage2Metabo/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo.svg)](https://github.com/AuReMe/metage2metabo/blob/main/LICENSE) [![Actions Status](https://github.com/AuReMe/metage2metabo/actions/workflows/pythonpackage.yml/badge.svg)](https://github.com/AuReMe/metage2metabo/actions/workflows/pythonpackage.yml) [![Documentation Status](https://readthedocs.org/projects/metage2metabo/badge/?version=latest)](https://metage2metabo.readthedocs.io/en/latest/?badge=latest) [![](https://img.shields.io/badge/doi-10.7554/eLife.61968-blueviolet.svg)](https://doi.org/10.7554/eLife.61968)

# M2M - metage2metabo
Metage2metabo is a Python3 (Python >= 3.8, tested with 3.8 and 3.9) tool to perform graph-based metabolic analysis starting from annotated genomes (**reference genomes or metagenome-assembled genomes**). It uses *Pathway Tools* in a automatic and parallel way to **reconstruct metabolic networks** for a large number of genomes. The obtained metabolic networks are then **analyzed individually and collectively** in order to get the **added value of metabolic cooperation in microbiota over individual metabolism** and to **identify and screen interesting organisms** among all.


m2m can be used as a whole workflow (``` m2m workflow ```, ``` m2m metacom ```) or steps can be performed individually (``` m2m recon ``` , ``` m2m iscope ``` , ``` m2m cscope ```, ``` m2m addedvalue ```, ``` m2m mincom ```, ``` m2m seeds ```).

A dedicated pipeline for the analysis of minimal community solution and their visualisation in power graphs is available with `m2m analysis`. More details in the [documentation](https://metage2metabo.readthedocs.io/en/latest/m2m_analysis.html).

**If you use M2M, please cite:**

* Belcour* A, Frioux* C, Aite M, Bretaudeau A, Hildebrand F, Siegel A. (2020). Metage2Metabo, microbiota-scale metabolic complementarity for the identification of key species. eLife 2020;9:e61968 [https://doi.org/10.7554/eLife.61968](https://doi.org/10.7554/eLife.61968).

* Frioux, C., Fremy, E., Trottier, C., & Siegel, A. (2018). Scalable and exhaustive screening of metabolic functions carried out by microbial consortia. Bioinformatics, 34(17), i934–i943. [https://doi.org/10.1093/bioinformatics/bty588](https://doi.org/10.1093/bioinformatics/bty588).

**If you use m2m_analysis, please cite <ins>additionally</ins>:**

* Bourneuf, L., & Nicolas, J. (2017). FCA in a Logical Programming Setting for Visualization-Oriented Graph Compression. Formal Concept Analysis, 14th International Conference, ICFCA 2017, Rennes, France, June 13-16, 2017, Proceedings. 89–105. [https://doi.org/10.1007/978-3-319-59271-8_6](https://doi.org/10.1007/978-3-319-59271-8_6).

**If you use m2m recon, please cite the appropriate Pathway Tools paper, in addition to the first two references.**

For a summary of M2M and its applications, you can take a look at these [poster-slides](https://hal.inria.fr/hal-03151934/document), presented during the [JOBIM 2020 conference](https://jobim2020.sciencesconf.org/?forward-action=index&forward-controller=index&lang=en).

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
  - [M2M commands](#m2m-commands)
  - [Analysis of the minimal solutions](#analysis-of-the-minimal-solutions)
  - [Release Notes](#release-notes)
  - [Additional features](#additional-features)
  - [Citation](#citation)
  - [Article data](#article-data)
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

This project is licensed under the GNU Lesser General Public License - see the [LICENSE.md](https://github.com/AuReMe/metage2metabo/blob/main/LICENSE) file for details.

## Documentation

A more detailled documentation is available at: [https://metage2metabo.readthedocs.io](https://metage2metabo.readthedocs.io/en/latest/).

## Technologies

Python 3 (Python 3.8 and 3.9 are tested). M2M uses a certain number of Python dependencies. An example of all these dependencies working for Ubuntu 18.04 is available in [requirements.txt](https://github.com/AuReMe/metage2metabo/blob/main/requirements.txt).
They can be installed with:
````sh
pip install -r requirements.txt --no-cache-dir
````
In particular, m2m relies on:
* [mpwt](https://github.com/AuReMe/mpwt) to automatize metabolic network reconstruction with Pathway Tools
* [padmet](https://github.com/AuReMe/padmet) to manage metabolic networks
* [menetools](https://github.com/cfrioux/MeneTools) to analyze individual metabolic capabilities using logic programming. **Requires MeneTools > 3.4**
* [miscoto](https://github.com/cfrioux/miscoto) to analyze collective metabolic capabilities and select communities within microbiota using logic programming. **Requires MiSCoTo > 3.2**

Also, m2m_analysis relies on other packages:
* [networkx](https://github.com/networkx/networkx) to create graph from miscoto results
* [ete3](https://github.com/etetoolkit/ete) to add taxonomy information on the graph if you used mpwt taxon file
* [powergrasp](https://github.com/Aluriak/PowerGrASP) to compress networkx graph

## Requirements

* [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) version 23.0 or higher (free for [academic users](https://biocyc.org/download-bundle.shtml)) is **required for m2m workflow and m2m recon**.  Metage2Metabo uses mpwt for multiprocessing and mpwt is not usable on Windows. Therefore, the reconstruction step of Metage2Metabo **is not available on Windows**.
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
        On Linux and MacOS: ```export PATH=$PATH:/your/install/directory/pathway-tools```. 
        Consider adding Pathway Tools in ```$PATH``` permanently by running
        ````sh
        echo 'export PATH="$PATH:your/install/directory/pathway-tools:"' >> ~/.bashrc
        ````
        Then source the bashrc file:
        ````sh
        source ~/.bashrc

* [Oog Power Graph Command line tool](https://github.com/AuReMe/metage2metabo/tree/main/external_dependencies/Oog_CommandLineTool2012) to create a svg file from the compressed graph at the end of m2m_analysis. This tool is a jar file (``Oog.jar``) so Java is needed to use it.

## Installation

Developed and tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14) with Python3.8.

Continuous Integration using GitHub Actions with Python3.8 and Python3.9 on ubuntu-latest, macos-latest and windows-latest ([corresponding virtual environment](https://docs.github.com/en/free-pro-team@latest/actions/reference/specifications-for-github-hosted-runners#supported-runners-and-hardware-resources)).

### Installation with pip

```
pip install Metage2Metabo
```

### Availability on Docker and Singularity

Due to Pathway-Tools license, Docker or Singularity images are not available publicly.

But you can create these images by using the Dockerfile and Singularity recipes available inside the recipes folder.
With these files, you can create container with Pathway-Tools and m2m.

More informations in the [Docker and Singularity Documentation](https://metage2metabo.readthedocs.io/en/latest/install.html#installation-with-docker).

## M2M commands

M2M commands are listed in the [Commands Documentation](https://metage2metabo.readthedocs.io/en/latest/command.html).

````
Copyright (C) Dyliss & Pleiade
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

## Analysis of the minimal solutions

M2M performs a community minimization to find the union and intersection of the minimal communities. But it is possible to analyze all the minimal communities.
M2M has a second command-line, named m2m_analysis that performs this analysis. This method is slower than m2m as all sollutions are enumerated.
Then it creates a solutions graph and compresses it in a powergraph. Then it creates visualization (html file and optionnaly svg files).

More information about this command in the [m2m_analysis Documentation](https://metage2metabo.readthedocs.io/en/latest/m2m_analysis.html).

````
usage: m2m_analysis [-h] [-v] {enum,graph,powergraph,workflow} ...

Detection of key species among communities.
 For specific help on each subcommand use: m2m_analysis {cmd} --help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid subcommands:

  {enum,graph,powergraph,workflow}
    enum                enumeration using miscoto
    graph               graph creation with enumeration solution
    powergraph          powergraph creation and visualization
    workflow            whole workflow

Optional: Oog.jar file (https://github.com/AuReMe/metage2metabo/tree/main/external_dependencies) for powergraph svg creation.
````

## Release Notes

Changes between version are listed on the [release page](https://github.com/AuReMe/metage2metabo/releases).

## Additional features

M2M relies on packages that can also be used independantly with more features:
* [mpwt](https://github.com/AuReMe/mpwt): command-line and multi-process solutions to run Pathway Tools. Suitable to multiple reconstruction, for example genomes of a microbiota
* [menetools](https://github.com/cfrioux/MeneTools): individual metabolic capabilities analysis using graph-based producibility criteria
* [miscoto](https://github.com/cfrioux/miscoto): community selection and metabolic screening in large-scal microbiotas, with or without taking a host into account

## Citation

If you use Metage2Metabo, please cite:

**Belcour\* A, Frioux\* C, Aite M, Bretaudeau A, Hildebrand F, Siegel A. Metage2Metabo, microbiota-scale metabolic complementarity for the identification of key species. eLife 2020;9:e61968 [https://doi.org/10.7554/eLife.61968](https://doi.org/10.7554/eLife.61968).**

Also when using ``m2m``, please cite the following articles:

- ``MeneTools`` for individual scope computation:

Aite M, Chevallier M, Frioux C, Trottier C, Got J, Cortés M P, Mendoza S N, Carrier G, Dameron O, Guillaudeux N, Latorre M, Loira N, Markov G V, Maass A, Siegel A. Traceability, reproducibility and wiki-exploration for “à-la-carte” reconstructions of genome-scale metabolic models. PLOS Computational Biology 2018;14:e1006146. [https://doi.org/10.1371/journal.pcbi.1006146](https://doi.org/10.1371/journal.pcbi.1006146).

- ``MiSCoTo`` for community scope computation and minimal community selection:

Frioux C, Fremy E, Trottier C, Siegel A. Scalable and exhaustive screening of metabolic functions carried out by microbial consortia. Bioinformatics 2018;34:i934–43. [https://doi.org/10.1093/bioinformatics/bty588](https://doi.org/10.1093/bioinformatics/bty588).

If you use ``m2m recon``, please cite additionally:

- ``Pathway Tools`` for the reconstruction of draft metabolic networks (the article can be not up-to-date, look at the [Publications](https://biocyc.org/publications.shtml) on the BioCyc site):

Karp P D, Midford P E, Billington R, Kothari A, Krummenacker M, Latendresse M, Ong W K, Subhraveti P, Caspi R, Fulcher C, Keseler I M, Paley SM. Pathway Tools version 23.0 update: software for pathway/genome informatics and systems biology. Briefings in Bioinformatics 2021;22:109–126. [https://doi.org/10.1093/bib/bbz104](https://doi.org/10.1093/bib/bbz104).

- ``padmet`` library for metabolic network storage (same article for ``MeneTools``):

Aite M, Chevallier M, Frioux C, Trottier C, Got J, Cortés M P, Mendoza S N, Carrier G, Dameron O, Guillaudeux N, Latorre M, Loira N, Markov G V, Maass A, Siegel A. Traceability, reproducibility and wiki-exploration for “à-la-carte” reconstructions of genome-scale metabolic models. PLOS Computational Biology 2018;14:e1006146. [https://doi.org/10.1371/journal.pcbi.1006146](https://doi.org/10.1371/journal.pcbi.1006146).


If you use ``m2m_analysis``, please cite additionally:

- ``networkx`` for graph solution creation:

Hagberg A A, Schult D A, Swart P J. Exploring Network Structure, Dynamics, and Function using NetworkX, in: Varoquaux, G., Vaught, T., Millman, J. (Eds.), . Presented at the Proceedings of the Python in Science Conference (SciPy) 2008. 11–15. [http://conference.scipy.org/proceedings/SciPy2008/paper_2/](http://conference.scipy.org/proceedings/SciPy2008/paper_2/)

- ``Powergrasp`` for power graph compression:

Bourneuf L, Nicolas J. FCA in a Logical Programming Setting for Visualization-Oriented Graph Compression. In: Bertet K, Borchmann D, Cellier P, Ferre´ S (Eds). ICFCA 2017: Formal Concept Analysis 2017. Springer. 89–105. [https://doi.org/10.1007/978-3-319-59271-8_6](https://doi.org/10.1007/978-3-319-59271-8_6).


- ``Oog command line tool`` for power graph visualisation:

Royer L, Reimann M, Andreopoulos B, Schroeder M, Unraveling Protein Networks with Power Graph Analysis. PLOS Computational Biology 2008;4:e1000108. [https://doi.org/10.1371/journal.pcbi.1000108](https://doi.org/10.1371/journal.pcbi.1000108).

- ``ete3`` for taxonomic information used in power graphs:

Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Molecular Biology and Evolution 2016;33:1635–1638. [https://doi.org/10.1093/molbev/msw046](https://doi.org/10.1093/molbev/msw046).


## Article data

Data used to create figures and tables are listed in the [article_data](https://github.com/AuReMe/metage2metabo/tree/main/article_data) folder, it contains:

- [gsmn_characteristics](https://github.com/AuReMe/metage2metabo/tree/main/article_data/gsmn_characteristics): scripts and tables to show the characteristics of draft metabolic networks created by M2M for gut, rumen and diabetes dataset.
- [diabetes_study](https://github.com/AuReMe/metage2metabo/tree/main/article_data/diabetes_study): scripts and tables to create the figures of the diabetes analyses in the article.

## Authors
[Clémence Frioux](https://cfrioux.github.io/) and [Arnaud Belcour](https://arnaudbelcour.github.io/blog/), `Univ Bordeaux, Inria, INRAE, Bordeaux, France`, `Univ Grenoble Alpes, Inria, Grenoble, France` and `Univ Rennes, Inria, CNRS, IRISA, Rennes, France`.

## Acknowledgement
People of Pathway Tools (SRI International) for their help integrating Pathway Tools with command line and multiprocessing in the [mpwt](https://github.com/AuReMe/mpwt) package, used in M2M.
