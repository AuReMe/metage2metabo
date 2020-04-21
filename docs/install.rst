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

        Install though a .dmg installer with a graphical interface.

        * **Warning**

        /!\\ For all OS, Pathway-Tools must be in ``$PATH``.

        On Linux and MacOS: ``export PATH=$PATH:your/install/directory/pathway-tools``.

        Consider adding Pathway Tools in ``$PATH`` permanently by running

        .. code:: sh

            echo 'export PATH="$PATH:your/install/directory/pathway-tools:"' >> ~/.bashrc

m2m relies on several Python packages:

* `mpwt <https://github.com/AuReMe/mpwt>`__ to automatize metabolic network reconstruction with Pathway-Tools.

* `padmet <https://github.com/AuReMe/padmet>`__ to manage metabolic networks.

* `menetools <https://github.com/cfrioux/MeneTools>`__ to analyze individual metabolic capabilities using logic programming.

* `miscoto <https://github.com/cfrioux/miscoto>`__ to analyze collective metabolic capabilities and select communities within microbiota using logic programming.

A list of Python dependencies is available in `requirements.txt <https://github.com/AuReMe/metage2metabo/blob/master/requirements.txt>`__

Tested on Linux (Ubuntu, Fedora, Debian) and MacOs (version 10.14).

Tested with Python v3.6, v3.7, v3.8.

Installation with pip
=====================

.. code:: sh

    pip install Metage2Metabo

Installation with Docker
========================

To create the m2m image, use the Dockerfile found in `Recipes <https://github.com/AuReMe/metage2metabo/tree/master/recipes>`__ of the Github repository. Note that the **Pathway-Tools installer** and the **Oog.jar file** need to be placed in the same folder than the Dockerfile.
The name of the installer file is currently hardcoded in the Dockerfile. Hence it must be changed if you use a different version of Pathway-Tools. Please note that the following commands (especially due to the use of root privileges) apply to Linux OS.

.. code:: sh

    # Launch docker.
    sudo systemctl start docker

    " Build image locally.
    sudo docker build -t my_image .

To create and launch the container in interactive mode:

.. code:: sh

    sudo docker run -ti -v /my/path/to/my/data:/shared --name="my_container" my_image bash

Then you can exit the contanier with ``exit``. You can launch again the container with:

.. code:: sh

    sudo docker start my_container

    sudo docker exec -ti my_container bash

Installation with Singularity (e.g. on a cluster)
=================================================

Singularity [Ku2017]_ can be used to launch m2m on a cluster. Please refer to the `recipe <https://github.com/AuReMe/metage2metabo/tree/master/recipes>`__   of the Github repository of the project.
The Singularity container has to be created from the recipe. You might need to do it on a personal computer since it requires administrator priviledges.
To use the container on a cluster, the path to Pathway Tools ptools folder should be indicated in the recipe. Therefore, you have to replace '/external/folder/ptools' with the path where you want to put the ptools-local folder (which will contain the PGDB created by Pathway-Tools).

Like for the Dockerfile, Pathway-Tools installer is hardcoded in the recipe so if you use another version, you have to modify the recipe.
And the **Pathway Tools installer** and **Oog.jar file** must be stored in the same folder than the Singularity recipe.

To create a container named m2m.sif:

.. code:: sh

    sudo singularity build m2m.sif Singularity

To use Pathway-Tools, a .ncbirc file is required in the home directory, containing the path to Blast:

.. code:: sh

    .ncbirc:

    [ncbi]\nData=/usr/bin/data

*Dealing with Pathway Tools ptools local folder*.
You might need an external ptools-local folder when working on a cluster. A solution is to create the ptools-local in a local folder then move it inside the Singularity container.
Eventually, you have to move it outside the Singularity container after it has been built.

First, enter the Singularity container and mount the external folder:

.. code:: sh

    singularity run -B /external/folder:/external/folder m2m.sif


Then move the ptools-local folder from the Singularity folder to the folder in your local environment.

.. code:: sh

    cp -r /opt/ptools-local /external/folder

This will move the ptools-local folder (with permissions) from Singularity container to the local machine.

In this way, PGDBs can be stored in the folder outside your container.

Finally, you can launch jobs with the Singularity container by giving a sh file containg m2m commands.

.. code:: sh

    m2m.sh:

    m2m workflow -g genomes_dir -s seeds.sbml -o output_dir -c cpu_number

So you can encapsulate it in a sh script:

.. code:: sh

    my_script.sh:

    #!/bin/bash

    # Don't forget to source the Singularity environment if needed.
    . /local/env/envsingularity.sh

    singularity exec m2m.sif bash m2m.sh

This file can now be launched on a cluster, for example with SLURM [Yo2003]_:

.. code:: sh

    sbatch --cpus-per-task=4 --mem=8G my_script.sh



.. [Ku2017] Kurtzer GM, Sochat V, Bauer MW (2017) Singularity: Scientific containers for mobility of compute. PLOS ONE 12(5): e0177459.

.. [Yo2003] Yoo, Andy B., Jette, Morris A., Grondona, Mark (2003). SLURM: Simple Linux Utility for Resource Management. Job Scheduling Strategies for Parallel Processing. Lecture Notes in Computer Science. 2862. p. 44. 
