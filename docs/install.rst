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

m2m relies on several Python packages:

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

To create the image using the Dockerfile in the github repository (in recipes), you need a Pathway-Tools installer placed in the same folder than the dockerfile.
The name of the file is currently hardcoded (so if you use another version it must be changed).


.. code:: sh

    # Launch docker.
    sudo systemctl start docker

    sudo docker build -t my_image .

To launch the container in interactive mode:

.. code:: sh

    sudo docker run -ti -v /my/path/to/my/data:/shared --name="my_container" my_image bash

Installation with Singularity
=============================

To launch m2m in a cluster we use Singularity.
First, we need to create the Singularity image locally.
But to use the image on the cluster, we have to put the cluster path in the Singularity recipe.
So you have to replace '/external/folder/ptools' with the path where you want to put the ptools-local folder (which contains PGDB for Pathway-Tools).

To create the image named m2m.sif (using the Singularity file in the recipes folder):

.. code:: sh

    sudo singularity build m2m.sif Singularity

To use Pathway-Tools, you need a file named .ncbirc in your home and containing the path to Blast:

.. code:: sh

    .ncbirc:

    [ncbi]\nData=/usr/bin/data

So in a cluster you need to create this file in your home.

To have an external ptools-local folder (mandatory when using the image on cluster), we have implemented an ugly hack.
The idea is that it creates the ptools-local in a local folder then it moves it inside the Singularity image.
So you have to move it outside the Singularity image after it has been built.

First, enter the Singularity image:

.. code:: sh

    singularity run m2m.sif


Then move the ptools-local folder from the Singularity folder to the folder in your local environment.

.. code:: sh

    cp -r /opt/ptools /external/folder/ptools

This will move the ptools-local folder (with permissions) from Singularity container to the local machine.

In this way, PGDBs can be stored in the folder outside your container.

Then you can launch jobs with the Singularity image by giving a sh file containg m2m commands.

.. code:: sh

    m2m.sh:

    m2m workflow -g genomes_dir -s seeds.sbml -o output_dir -c cpu_number

So you can encapsulate it in a sh script:

.. code:: sh

    my_script.sh:

    #!/bin/bash

    # Don't forget to source the Singularity environment.
    . /local/env/envsingularity.sh

    singularity exec m2m.sif bash m2m.sh

This file can now be launched on a cluster, for example (in SLURM):

.. code:: sh

    sbatch --cpus-per-task=4 --mem=8G my_script.sh
