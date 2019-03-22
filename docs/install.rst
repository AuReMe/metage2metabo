=============================
Requirements and Installation
=============================

Requirements
============

* `Pathway-Tools <http://bioinformatics.ai.sri.com/ptools/>`__ version 22.5 or higher (free for `academic users <https://biocyc.org/download-bundle.shtml>`__)
    * Pathway-Tools requirements
        * **Linux**: Gnome terminal and Libxm4

        .. code:: sh

            apt-get update && apt-get install gnome-terminal libxm4

        * **All OS**: `NCBI Blast <https://www.ncbi.nlm.nih.gov/books/NBK279671/>`__ and a ncbirc file in user's home directory
            * Install with apt-get

            .. code:: sh

                apt-get update && apt-get install gnome-terminal libxm4 ncbi-blast+
                echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

            * Install with a dmg installer on MacOS
    * Pathway-Tools install
        * **Linux**

        .. code:: sh

            chmod +x ./pathway-tools-22.5-linux-64-tier1-install
            ./pathway-tools-22.5-linux-64-tier1-install

        and follow the instructions during the interactive install

        *For a silent install*:

        .. code:: sh

            ./pathway-tools-22.5-linux-64-tier1-install --InstallDir your/install/directory/pathway-tools --PTOOLS_LOCAL_PATH your/chosen/directory/for/data/ptools --InstallDesktopShortcuts 0 --mode unattended

        * **MacOS**

        Dmg installer with a graphical interface.

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).
Tested with Python3.6

Installation with pip
=====================

.. code:: sh

    pip install m2m

Installation with docker?
=========================

Installation with Singularity
=============================