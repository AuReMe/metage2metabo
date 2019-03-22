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

            * Install with a `dmg installer <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`__ on MacOS
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

        * **Warning**

        /!\\ For all OS, Pathway-Tools must be in ``$PATH``.

        On Linux and MacOS: ``export PATH=$PATH:your/install/directory/pathway-tools``.

        Consider adding Pathway Tools in ``$PATH`` permanently by running

        .. code:: sh

            echo 'export PATH="$PATH:your/install/directory/pathway-tools:"' >> ~/.bashrc

m2m relies on python package:

* `mpwt <https://github.com/AuReMe/mpwt>`__ to automatize metabolic network reconstruction with Pathway-Tools.

* `padmet <https://github.com/AuReMe/padmet>`__ and `padmet-utils <https://github.com/AuReMe/padmet-utils>`__ to manage metabolic networks.

* `menetools <https://github.com/cfrioux/MeneTools>`__ to analyze individual metabolic capabilities using logic programming.

* `miscoto <https://github.com/cfrioux/miscoto>`__ to analyze collective metabolic capabilities and select communities within microbiota using logic programming.

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).

Tested with Python3.6

Installation with pip
=====================

.. code:: sh

    pip install m2m

Installation with Docker
========================

Installation with Singularity
=============================
