# M2M - metage2metabo

Metage2metabo is a Python3 workflow to perform graph-based metabolic analysis starting from annotated genomes. It uses Pathway Tools in a automatic and parallel way to reconstruct metabolic networks for a large number of genomes. The obtained metabolic networks are then analyzed individually and collectively in order to get the added value of cooperation in microbiota over individual metabolism, and to identify and screen interesting organisms among all

## Install

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).
Tested with Python3.6

### Install with pip

```pip install metage2metabo```

### Availability on Docker and Singularity

```details ```

## Requirements

* [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) version 22.5 or higher (free for [academic users](https://biocyc.org/download-bundle.shtml))
    * Pathway Tools requirements
        * **Linux**: Gnome terminal and Libxm4
        ```
        apt-get update && apt-get install gnome-terminal libxm4
         ```
        * **All OS**: [NCBI Blast](https://www.ncbi.nlm.nih.gov/books/NBK279671/) and a ncbirc file in user's home directory
            * Install with apt-get
            ```sh
            apt-get update && apt-get install gnome-terminal libxm4 ncbi-blast+ 
            echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc
            ```
            * Install with a dmg installer on MacOS
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
