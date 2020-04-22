============
m2m Tutorial
============
Test data is avaible in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/test>`__.
It contains enough data to run the different subcommands.

Outputs and logs
-----------------

By default, ``m2m`` writes into the console (stdout) and into a log file located at the root of the results directory and named after the subcommand that wad executed. The option ``-q`` can be given to any ``m2m`` subcommand to write in the log file and write to stdout only the warnings, errors and critical issues, together with the path to the log file.

m2m recon
---------
``m2m recon`` runs metabolic network reconstruction for all annotated genomes, using Pathway Tools. It can be done with multiple CPUs, in which case the number of allocated/available CPU has to be given as a optional argument.

It uses the following mandatory inputs (run ``m2m recon --help`` for optional arguments):

-g directory           directory of annotated genomes
-o directory           output directory for results

Optional arguments:

-c int           number of CPU for multi-processing
--clean          option to rerun every reconstruction 
                 even if found in ptools-local
--noorphan       ignore the reactions without gene or 
                 protein association in final metabolic networks
-p               create padmet files from PGDB
-l int           specify the level for the sbml to be created
-q               quiet mode

The input genomic data has to follow a strict structure:

::

    input_folder
        ├── organism_1
        │   └── organism_1.gbk
        ├── organism_2
        │   └── organism_2.gbk
        ├── organism_3
        │   └── organism_3.gbk
        ..
        └── organism_n         
            └── organism_n.gbk

This structure can be observed in the `workflow_genomes.tar.gz` file in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/metage2metabo/workflow_data>`__.

By extracting this file, you will find the

::

    workflow_genomes
        ├── GCA_003433665
        │   └── GCA_003433665.gbk
        ├── GCA_003433675
        │   └── GCA_003433675.gbk

In addition, the genbank files also follow the strict requirements of Pathway Tools.
Please go check the `documentation of mpwt <https://github.com/AuReMe/mpwt#genbank>`__ for more details, especially if you get errors from Pathway Tools at the ``recon`` step.

.. code:: sh

    m2m recon -g workflow_genomes -o output_directory -c cpu_number [--clean] [--orphan] [-p] [-l sbml_level]

* standard output
    .. code:: 

        ######### Running metabolic network reconstruction with Pathway Tools #########
        ~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~
        Checking inputs for GCA_003433675: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for GCA_003433675.
        Checking inputs for GCA_003433665: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for GCA_003433665.
        ----------End of creation of input data from Genbank/GFF/PF: 0.18s----------
        ~~~~~~~~~~Inference on the data~~~~~~~~~~
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho test/GCA_003433675/
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho test/GCA_003433665/
        ~~~~~~~~~~Check inference~~~~~~~~~~

        2 builds have passed!

        ----------End of PathoLogic inference: 323.64s----------
        ~~~~~~~~~~Creation of the .dat files~~~~~~~~~~
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load test/GCA_003433665/dat_creation.lisp
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load test/GCA_003433675/dat_creation.lisp
        ~~~~~~~~~~Check .dat~~~~~~~~~~
        gca_003433665cyc: 23 out of 23 dat files create.
        gca_003433675cyc: 23 out of 23 dat files create.
        ----------End of dat files creation: 158.31s----------
        ~~~~~~~~~~End of Pathway Tools~~~~~~~~~~
        ~~~~~~~~~~Moving result files~~~~~~~~~~
        ----------End of moving fimes: 0.12s----------
        ----------mpwt has finished in 482.34s! Thank you for using it.
        ######### Creating SBML files #########
        ######### Stats GSMN reconstruction #########
        Number of genomes: 2
        Number of reactions in all GSMN: 2299
        Number of compounds in all GSMN: 2446
        Average reactions per GSMN: 1654.00(+/- 743.88)
        Average compounds per GSMN: 1840.50(+/- 723.37)
        Average genes per GSMN: 831.00(+/- 435.58)
        Percentage of reactions associated with genes: 78.17(+/- 3.23)
        --- Recon runtime 486.88 seconds ---

        PGDB created in out1_test/pgdb
        SBML files created in out1_test/sbml
        --- Total runtime 486.89 seconds ---

        The output shows that PGDB are created with Pathway Tools. Then the .dat files are extracted and used to build SBML files of the metabolic models.
* files outputs
    * In `output_directory/pgdb`, the .dat files of Pathway Tools. The corresponding SBMLs are in `output_directory/sbml`. The structure of the output directory after this ``recon`` command is shown below :

    ::

        output_directory/
        ├── m2m_recon.log
        ├── pgdb
        │   ├── GCA_003433665
        │   │   ├── classes.dat
        │   │   ├── compound-links.dat
        │   │   ├── compounds.dat
        │   │   ├── dnabindsites.dat
        │   │   ├── enzrxns.dat
        │   │   ├── gene-links.dat
        │   │   ├── genes.dat
        │   │   ├── pathway-links.dat
        │   │   ├── pathways.dat
        │   │   ├── promoters.dat
        │   │   ├── protein-features.dat
        │   │   ├── protein-links.dat
        │   │   ├── proteins.dat
        │   │   ├── protligandcplxes.dat
        │   │   ├── pubs.dat
        │   │   ├── reaction-links.dat
        │   │   ├── reactions.dat
        │   │   ├── regulation.dat
        │   │   ├── regulons.dat
        │   │   ├── rnas.dat
        │   │   ├── species.dat
        │   │   ├── terminators.dat
        │   │   └── transunits.dat
        │   └── GCA_003433675
        │       ├── classes.dat
        │       ├── compound-links.dat
        │       ├── compounds.dat
        │       ├── dnabindsites.dat
        │       ├── enzrxns.dat
        │       ├── gene-links.dat
        │       ├── genes.dat
        │       ├── pathway-links.dat
        │       ├── pathways.dat
        │       ├── promoters.dat
        │       ├── protein-features.dat
        │       ├── protein-links.dat
        │       ├── proteins.dat
        │       ├── protligandcplxes.dat
        │       ├── pubs.dat
        │       ├── reaction-links.dat
        │       ├── reactions.dat
        │       ├── regulation.dat
        │       ├── regulons.dat
        │       ├── rnas.dat
        │       ├── species.dat
        │       ├── terminators.dat
        │       └── transunits.dat
        └── recon_stats.tsv
        └── sbml
            ├── GCA_003433665.sbml
            └── GCA_003433675.sbml

        * Finally, in the input directory, some files are also generated automatically by Pathway Tools
        ::
            
            recon_data/
            ├── GCA_003433665
            │   ├── dat_creation.lisp
            │   ├── GCA_003433665.gbk
            │   ├── genetic-elements.dat
            │   ├── organism-params.dat
            │   └── pathologic.log
            └── GCA_003433675
                ├── dat_creation.lisp
                └── GCA_003433675.gbk
                ├── genetic-elements.dat
                ├── organism-params.dat
                ├── pathologic.log


m2m iscope, cscope and addedvalue
---------------------------------
The three subcommands require metabolic networks under the SBML format. Some metabolic networks are available as a compressed archive in `metabolic_data`. Uncompress the file and the directory can be fed to the subcommands. These commands also require a seeds file comprising the metabolic compounds available to assess reachability/producibility in the models. This seeds file needs to be in SBML format. You can use the one in the `metabolic data` directory.

Optional: create the seeds SBML file
*************************************
To create a seeds file starting from a list of metabolic identifiers (matching identifiers of compounds of the organisms metabolic networks), you can use the ``m2m seeds`` command:

.. code:: sh

    m2m seeds --metabolites metabolites_file.txt -o output/directory

The resulting seeds file will be created in output/directory/seeds.sbml

An example of structure of the metabolites file is the following:

.. code:: 

    M_AMMONIA_c
    M_ZN__43__2_c
    M_CARBON__45__DIOXIDE_c
    M_OXYGEN__45__MOLECULE_c

The resulting SBML will have such a design:

.. code:: xml

    <?xml version="1.0" encoding="UTF-8"?>
        <sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1">
        <model id="metabolites">
            <listOfSpecies>
            <species id="M_AMMONIA_c" name="AMMONIA" compartment="c"/>
            <species id="M_ZN__43__2_c" name="ZN+2" compartment="c"/>
            <species id="M_CARBON__45__DIOXIDE_c" name="CARBON-DIOXIDE" compartment="c"/>
            <species id="M_OXYGEN__45__MOLECULE_c" name="OXYGEN-MOLECULE" compartment="c"/>
            </listOfSpecies>
    </model>
    </sbml>

iscope
*******

It uses the following mandatory inputs (run ``m2m iscope --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional argument
-q                     quiet mode

.. code:: sh

    m2m iscope -n toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory/

* standard output
    .. code:: 

        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json
        17 metabolic models considered.
        135 metabolites in core reachable by all organisms (intersection)
        625 metabolites reachable by individual organisms altogether (union), among which 93 seeds (growth medium)
        max metabolites in scope 477
        min metabolites in scope 195
        average number of metabolites in scope 308.71 (±82.59)

    These results mean that 135 metabolites can be reached by all organisms. When gathering reachable metabolites for all organisms, the union consists of 625 metabolites (including the seeds). Finally metrics show the min, max and average number of compounds in all scopes
* files outputs
    * In `output_directory/indiv_scopes/indiv_scopes.json`. A json file that can be easily loaded as a dictionary (or humanly read as it it) that contains the set of reachable metabolites for each organism. /!\\ Warning: the seeds are included in the scopes, hence they will never be empty. Logs are written in `output_directory/m2m_iscope.log` .

cscope
*******

It uses the following mandatory inputs (run ``m2m cscope --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file                targets SBML file
-o directory           output directory for results
-m file                host metabolic network SBML file

Optional arguments:

-m file                host metabolic network SBML file
-q                     quiet mode

.. code:: sh

    m2m cscope -n toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory/

* standard output
    .. code::

        ######### Creating metabolic instance for the whole community #########
        Created instance in output_directory/community_analysis/miscoto_om6hubmz.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        651 metabolites reachable by the whole community/microbiota:
        M_CPD__45__5802_c, M_XANTHOSINE__45__5__45__PHOSPHATE_c, M_INDOLEYL__45__CPD_c, M_CPD__45__470_c, M_5__45__HYDROXYISOURATE_c, [...]

    651 metabolites are reachable by the microbiota. This does not include the seeds. The list of metabolites is given in output. 
* files outputs
    * In addition to the logs at the root of the results directory, a json file with the results is created in `output_directory/community_analysis/indiv_scopes.json`.

addedvalue
**********

``m2m addedvalue`` uses the previously two subcommands to compute the added value of combining metabolisms in the microbiota (i.e. consider metabolic cooperation) with respect to studying individually the metabolism of each organism. 
It uses the following mandatory inputs (run ``m2m addedvalue --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional arguments:

-m file                host metabolic network SBML file
-q                     quiet mode

.. code:: sh

    m2m addedvalue -n toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory/

* standard output
    .. code::

        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json
        17 metabolic models considered.
        135 metabolites in core reachable by all organisms (intersection)
        625 metabolites reachable by individual organisms altogether (union), among which 93 seeds (growth medium)
        max metabolites in scope 477
        min metabolites in scope 195
        average number of metabolites in scope 308.71 (±82.59)
        M_D__45__RIBULOSE__45__1__45__P_c, M_ISOGLUTAMINE_c, M_RIBULOSE__45__5P_c, M_MET_c, M_CPD__45__10775_c, M_DGDP_c, M_5__45__PHOSPHO__45__RIBOSYL__45__GLYCINEAMIDE_c, M_ADENYLOSUCC_c, M_ISOCHORISMATE_c, [...]
        ######### Creating metabolic instance for the whole community #########
        Created instance in output_directory/community_analysis/miscoto_j9khdvzz.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        651 metabolites reachable by the whole community/microbiota:
        M_D__45__RIBULOSE__45__1__45__P_c, M_ISOGLUTAMINE_c, M_RIBULOSE__45__5P_c, M_CPD__45__10775_c, M_DGDP_c, M_5__45__PHOSPHO__45__RIBOSYL__45__GLYCINEAMIDE_c, M_OH__45__HEXANOYL__45__COA_c, M_ADENYLOSUCC_c,[...]
        Added value of cooperation over individual metabolism: 119 newly reachable metabolites:
        M_OH__45__HEXANOYL__45__COA_c, M_CPD__45__12307_c, M_CPD__45__12173_c, M_2__45__METHYL__45__ACETO__45__ACETYL__45__COA_c, [...]
        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml

    As you can see, the individual and community scopes are run again. In addition to the previous outputs, the union of all individual scopes and the community scopes are printed. Finally, the difference between the two sets, that is to say the metabolites that can only be produced collectively (i.e. by at least two bacteria cooperating) is displayed. Here it consists of 119 metabolites. 
* files outputs
    * A targets SBML file is generated. It can be used with `` m2m mincom`` . Newly producible metabolites are written in a json file. The json files associated to ``iscope`` and ``cscope`` are also produced.

    ::

        output_directory/
        ├── m2m_addedvalue.log
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── addedvalue.json
        │   └── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json


m2m mincom
----------
`m2m mincom` requires an additional target file that is available in `metabolic_data` or can be generated by `m2m addedvalue` in which case it will be stored in `result_directory/community_analysis/targets.sbml`

It uses the following mandatory inputs (run ``m2m mincom --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file                targets SBML file
-o directory           output directory for results

Optional arguments:

-m file                host metabolic network SBML file
-q                     quiet mode

.. code:: sh

    m2m mincom -n toy_bact -s metabolic_data/seeds_toy.sbml -t metabolic_data/targets_toy.sbml -o output_directory/

* standard output
    .. code::

        ######### Creating metabolic instance for the whole community #########
        Created instance in output_directory/community_analysis/miscoto_36t8lqe_.lp
        Running minimal community selection
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities = 13
        GCA_003437375
        GCA_003437945
        GCA_003437195
        GCA_003437295
        GCA_003437815
        GCA_003437595
        GCA_003437885
        GCA_003437905
        GCA_003437715
        GCA_003437255
        GCA_003437055
        GCA_003437665
        GCA_003438055
        ######### Keystone species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Keystone species = 17
        GCA_003437195
        GCA_003437175
        GCA_003437945
        GCA_003437785
        GCA_003437295
        GCA_003437885
        GCA_003437715
        GCA_003437345
        GCA_003437255
        GCA_003437375
        GCA_003437325
        GCA_003437815
        GCA_003437595
        GCA_003437905
        GCA_003437055
        GCA_003437665
        GCA_003438055
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal community enabling the producibility of the target metabolites given as inputs
        Essential symbionts = 12
        GCA_003437375
        GCA_003437195
        GCA_003437295
        GCA_003437815
        GCA_003437595
        GCA_003437885
        GCA_003437905
        GCA_003437715
        GCA_003437255
        GCA_003437055
        GCA_003437665
        GCA_003438055
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal community enabling the producibility of the target metabolites given as inputs
        Alternative symbionts = 5
        GCA_003437325
        GCA_003437345
        GCA_003437175
        GCA_003437945
        GCA_003437785


    This output gives the result of minimal community selection. It means that for producing the 119 metabolic targets, a minimum of 13 bacteria out of the 17 is required. One example of such minimal community is given. In addition, the whole space of solution is studied. All bacteria (17) occur in at least one minimal community (keystone species). Finally, the intersection gives the following information: a set of 12 bacteria occurs in each minimal communtity. This means that these 12 bacteria are needed in any case (essential symbionts), and that any of the remaining 5 bacteria (alternative symbionts) can complete the missing function(s).
* files outputs
    * As for other commands, a json file with the results is produced in ``output_directory/community_analysis/comm_scopes.json``, together with logs at the root of the results directory.

m2m metacom
------------
`m2m metacom` runs all analyses: individual scopes, community scopes, and minimal community selection based on the metabolic added-value of the microbiota.

It uses the following mandatory inputs (run ``m2m metacom --help`` for optional arguments):

-n directory           directory of metabolic networks,
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional arguments:

-m file                host metabolic network SBML file
-q                     quiet mode

.. code:: sh

    m2m metacom -n metabolic_data/toy_bact -s metabolic_data/seeds_toy.sbml  -o output_directory

* standard output
    .. code::

        At least one SBML has not a suitable level for the tools. They will be transformed and created in output_directory/new_sbml/. The others will be copied in this directory
        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json
        17 metabolic models considered.

        135 metabolites in core reachable by all organisms (intersection)

        ...

        625 metabolites reachable by individual organisms altogether (union), among which 93 seeds (growth medium)

        ...

        intersection of scope 135
        union of scope 625
        max metabolites in scope 477
        min metabolites in scope 195
        average number of metabolites in scope 308.71 (+/- 82.59)
        --- Indiv scopes runtime 5.78 seconds ---

        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/metage2metabo/test/output_directory/community_analysis/miscoto_5iys6bfh.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        --- Community scope runtime 3.26 seconds ---


        Added value of cooperation over individual metabolism: 119 newly reachable metabolites:

        ...

        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml
        Setting these 119 as targets
        Running minimal community selection
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities = 13
        GCA_003437715
        GCA_003437665
        GCA_003437055
        GCA_003437375
        GCA_003437595
        GCA_003437195
        GCA_003437295
        GCA_003437255
        GCA_003437885
        GCA_003438055
        GCA_003437815
        GCA_003437905
        GCA_003437945
        ######### Keystone species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Keystone species = 17
        GCA_003437715
        GCA_003437665
        GCA_003437055
        GCA_003437375
        GCA_003437195
        GCA_003437295
        GCA_003437255
        GCA_003437785
        GCA_003438055
        GCA_003437325
        GCA_003437905
        GCA_003437945
        GCA_003437815
        GCA_003437595
        GCA_003437885
        GCA_003437345
        GCA_003437175
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal community enabling the producibility of the target metabolites given as inputs
        Essential symbionts = 12
        GCA_003437715
        GCA_003437665
        GCA_003437055
        GCA_003437375
        GCA_003437595
        GCA_003437195
        GCA_003437295
        GCA_003437255
        GCA_003437885
        GCA_003438055
        GCA_003437815
        GCA_003437905
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal community enabling the producibility of the target metabolites given as inputs
        Alternative symbionts = 5
        GCA_003437945
        GCA_003437785
        GCA_003437345
        GCA_003437175
        GCA_003437325
        --- Mincom runtime 2.28 seconds ---

        --- Total runtime 16.21 seconds ---

* files outputs
    * Files are created in the output_directory: the logs, json files with the results, targets in SBML.

    .. code ::

        output_directory/
        ├── m2m_metacom.log
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── mincom.json
        │   ├── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json

m2m workflow
------------
`m2m workflow` starts from metabolic network reconstruction and runs all analyses: individual scopes, community scopes, and minimal community selection based on the metabolic added-value of the microbiota.

It uses the following mandatory inputs (run ``m2m workflow --help`` for optional arguments):

-g directory           directory of annotated genomes
-s file                seeds SBML file
-o directory           output directory for results

Optional arguments:

-c int           number of CPU for multi-processing
--clean          option to rerun every reconstruction 
                 even if found in ptools-local
--noorphan       ignore the reactions without gene or 
                 protein association in final metabolic networks
-p               create padmet files from PGDB
-q               quiet mode

You can run the workflow analysis with the two genbanks files available in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/metage2metabo>`__ (`workflow_data`). Two genomes are available in the compressed archive workflow_genomes.tar.gz. The archive has to be uncompressed before testing.

.. code:: sh

    m2m workflow -g workflow_genomes -s workflow_data/seeds_workflow.sbml -o output_directory/

Or you can run the test argument (which use the same data):

Which uses the following mandatory inputs (run ``m2m test --help`` for optional arguments):

-o directory           output directory path

Optional arguments:

-q               quiet mode
-c int           cpu number for multi-processing

.. code:: sh

    m2m test -o output_directory

* standard outputs

    .. code ::

        ######### Running metabolic network reconstruction with Pathway Tools #########
        ~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~
        Checking inputs for GCA_003433675: missing dat_creation.lisp; genetic-elements.dat; organism-params.dat. Inputs file created for GCA_003433675.
        Checking inputs for GCA_003433665: missing dat_creation.lisp; genetic-elements.dat; organism-params.dat. Inputs file created for GCA_003433665.
        ----------End of creation of input data from Genbank/GFF/PF: 0.18s----------
        ~~~~~~~~~~Inference on the data~~~~~~~~~~
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho test//GCA_003433675/
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho test//GCA_003433665/
        ~~~~~~~~~~Check inference~~~~~~~~~~
        2 builds have passed!
        ----------End of PathoLogic inference: 367.75s----------
        ~~~~~~~~~~Creation of the .dat files~~~~~~~~~~
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load test//GCA_003433675/dat_creation.lisp
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load test//GCA_003433665/dat_creation.lisp
        ~~~~~~~~~~Check .dat~~~~~~~~~~
        gca_003433675cyc: 23 out of 23 dat files create.
        gca_003433665cyc: 23 out of 23 dat files create.
        ----------End of dat files creation: 162.97s----------
        ~~~~~~~~~~End of Pathway Tools~~~~~~~~~~
        ~~~~~~~~~~Moving result files~~~~~~~~~~
        ----------End of moving fimes: 0.19s----------
        ----------mpwt has finished in 531.10s! Thank you for using it.
        ######### Creating SBML files #########
        ######### Stats GSMN reconstruction #########
        Number of genomes: 2
        Number of reactions in all GSMN: 2026
        Number of compounds in all GSMN: 2095
        Average reactions per GSMN: 1437.00(+/- 678.82)
        Average compounds per GSMN: 1560.00(+/- 615.18)
        Average genes per GSMN: 893.00(+/- 475.18)
        Percentage of reactions associated with genes: 79.90(+/- 3.20)
        --- Recon runtime 535.64 seconds ---
        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in out/indiv_scopes/indiv_scopes.json
        2 metabolic models considered.
        123 metabolites in core reachable by all organisms (intersection)
        M_SULFATE_c M_DIMETHYL__45__D__45__RIBITYL__45__LUMAZINE_c M_CPD0__45__2472_c M_AMMONIUM_c M_MN__43__2_c M_CPD__45__10809_c M_7__45__CYANO__45__7__45__DEAZAGUANINE_c M_CPD__45__69_c M_H2CO3_c M_CPD__45__602_c M_CARBAMOYL__45__P_c M_NADP_c M_NADPH_c M_P3I_c M_L__45__RIBULOSE__45__5__45__P_c M_ADP_c M_PHOSPHORIBOSYL__45__ATP_c M_GUANINE_c M_CPD0__45__2474_c M_ALPHA__45__GLUCOSE_c M_GLC_c M_FE__43__3_c M_NA__43___c M_FE__43__2_c M_CPD__45__18238_c M_DIHYDRO__45__NEO__45__PTERIN_c M_CA__43__2_c M_GLYCOLLATE_c M_CPD__45__18085_c M_PHOSPHORIBULOSYL__45__FORMIMINO__45__AICAR__45__P_c M_FRUCTOSE__45__6P_c M_CPD0__45__1699_c M_AMP_c M_DPG_c M_GLYCEROL__45__3P_c M_7__45__AMINOMETHYL__45__7__45__DEAZAGUANINE_c M_GLC__45__1__45__P_c M_CPD__45__3_c M_AMINO__45__RIBOSYLAMINO__45__1H__45__3H__45__PYR__45__DIONE_c M_GUANOSINE__45__5DP__45__3DP_c M_DIHYDRONEOPTERIN__45__P3_c M_ATP_c M_RIBULOSE__45__5P_c M_DIHYDROXYACETONE_c M_GMP_c M_CPD__45__653_c M_ACETALD_c M_MG__43__2_c M_DGTP_c M_DIHYDROXY__45__BUTANONE__45__P_c M_NADH_c M_D__45__glucopyranose__45__6__45__phosphate_c M_PROTON_c M_FAD_c M_URATE_c M_CPD__45__13469_c M_DATP_c M_XANTHOSINE_c M_FORMATE_c M_CPD__45__15709_c M_XYLULOSE__45__5__45__PHOSPHATE_c M_Glucopyranose_c M_IMIDAZOLE__45__ACETOL__45__P_c M_CPD__45__14133_c M_Pi_c M_WATER_c M_FMN_c M_CELLOBIOSE_c M_CU__43___c M_CPD__45__15818_c M_INOSINE_c M_GDP__45__TP_c M_ZN__43__2_c M_GUANOSINE_c M_IMP_c M_DIHYDRONEOPTERIN__45__P_c M_HYPOXANTHINE_c M_ADENOSINE_c M_NAD_c M_RIBOSE__45__5P_c M_AICAR_c M_3__45__P__45__HYDROXYPYRUVATE_c M_RIBOSE__45__1P_c M_CPD__45__13043_c M_PHOSPHORIBOSYL__45__FORMIMINO__45__AICAR__45__P_c M_PROTON_e M_CO__43__2_c M_AMMONIA_c M_GLYCOLALDEHYDE_c M_G3P_c M_CPD0__45__1108_c M_CL__45___c M_DIAMINO__45__OH__45__PHOSPHORIBOSYLAMINO__45__PYR_c M_GDP_c M_GAP_c M_CPD__45__10330_c M_GTP_c M_PPI_c M_XANTHINE_c M_K__43___c M_FRUCTOSE__45__16__45__DIPHOSPHATE_c M_ADENINE_c M_CPD__45__1086_c M_DIHYDROXY__45__ACETONE__45__PHOSPHATE_c M_DIHYDROPTERIN__45__CH2OH__45__PP_c M_PRPP_c M_HCO3_c M_CU__43__2_c M_RIBOFLAVIN_c M_NITRATE_c M_PHOSPHORIBOSYL__45__AMP_c M_3OH__45__4P__45__OH__45__ALPHA__45__KETOBUTYRATE_c M_D__45__Ribofuranose_c M_XANTHOSINE__45__5__45__PHOSPHATE_c M_AMINO__45__OH__45__HYDROXYMETHYL__45__DIHYDROPTERIDINE_c M_CARBAMATE_c M_ERYTHRONATE__45__4P_c M_D__45__Ribopyranose_c M_ERYTHROSE__45__4P_c M_CO3_c M_D__45__SEDOHEPTULOSE__45__7__45__P_c M_CARBON__45__DIOXIDE_c M_D__45__ERYTHRO__45__IMIDAZOLE__45__GLYCEROL__45__P_c
        325 metabolites reachable by individual organisms altogether (union), among which 26 seeds (growth medium)
        M_APS_c M_CPD__45__11770_c M_ISOCHORISMATE_c M_PYRIDOXAL_c M_DIMETHYL__45__D__45__RIBITYL__45__LUMAZINE_c M_ETOH_c M_2__45__KETO__45__3__45__DEOXY__45__D__45__GLUCARATE_c M_AMMONIUM_c M_MAL_c M_CPD__45__10809_c M_7__45__CYANO__45__7__45__DEAZAGUANINE_c M_CPD__45__602_c M_NADP_c M_GLYOX_c M_4__45__IMIDAZOLONE__45__5__45__PROPIONATE_c M_IMINOASPARTATE_c M_ISOGLUTAMINE_c M_2__45__PG_c M_2__45__KETOGLUTARATE_c M_ADP_c M_CPD__45__9924_c M_ALPHA__45__GLUCOSE_c M_GLC_c M_UROCANATE_c M_CPD__45__13118_c M_FE__43__2_c M_CA__43__2_c M_ARABINOSE__45__5P_c M_GLYCOLLATE_c M_HYDROGEN__45__MOLECULE_c M_FORMAMIDE_c M_CPD__45__18085_c M_ADP__45__D__45__GLUCOSE_c M_AMP_c M_ENTEROBACTIN_c M_INDOLE_ACETATE_AUXIN_c M_ADP__45__L__45__GLYCERO__45__D__45__MANNO__45__HEPTOSE_c M_INDOLE_PYRUVATE_c M_GDP__45__4__45__DEHYDRO__45__6__45__DEOXY__45__D__45__MANNOSE_c M_PYRIDOXAL_PHOSPHATE_c M_CPD__45__4841_c M_4__45__PHOSPHONOOXY__45__THREONINE_c M_AMINO__45__RIBOSYLAMINO__45__1H__45__3H__45__PYR__45__DIONE_c M_PYRIDOXINE__45__5P_c M_CPD__45__14443_c M_L__45__ASPARTATE_c M_CPD__45__19753_c M_DIHYDROXYACETONE_c M_2__45__KETO__45__ISOVALERATE_c M_THREO__45__DS__45__ISO__45__CITRATE_c M_L__45__GLYCERALDEHYDE__45__3__45__PHOSPHATE_c M_PYRUVATE_c M_CPD__45__653_c M_ACETALD_c M_MG__43__2_c M_DIHYDROXY__45__BUTANONE__45__P_c M_CPD__45__13357_c M_NITRITE_c M_TARTRONATE__45__S__45__ALD_c M_SERYL__45__AMP_c M_NADH_c M_CPD0__45__2483_c M_CIT_c M_DEOXYGUANOSINE_c M_C__45__DI__45__GMP_c M_PYRIDOXINE_c M_CPD0__45__1905_c M_TYR_c M_4__45__hydroxybenzoate_c M_CPD__45__12367_c M_URATE_c M_CPD__45__13469_c M_DATP_c M_CPD__45__13851_c M_XANTHOSINE_c M_FORMATE_c M_1__45__AMINO__45__PROPAN__45__2__45__ONE__45__3__45__PHOSPHATE_c M_CPD__45__15709_c M_XYLULOSE__45__5__45__PHOSPHATE_c M_GLUCOSAMINE__45__1P_c M_IMIDAZOLE__45__ACETOL__45__P_c M_DEHYDROQUINATE_c M_CPD__45__14133_c M_WATER_c M_FMN_c M_CPD__45__13559_c M_CELLOBIOSE_c M_KDO__45__8P_c M_CU__43___c M_CPD__45__15818_c M_INOSINE_c M_CHORISMATE_c M_GUANOSINE_c M_ADENYLOSUCC_c M_IMP_c M_INDOLE_c M_NAD_c M_ZN__43__2_e M_RIBOSE__45__5P_c M_O__45__SUCCINYLBENZOATE_c M_GDP__45__4__45__DEHYDRO__45__6__45__L__45__DEOXYGALACTOSE_c M_MANNOSE__45__1P_c M_DEOXY__45__RIBOSE__45__5P_c M_DEOXY__45__D__45__RIBOSE__45__1__45__PHOSPHATE_c M_SUPER__45__OXIDE_c M_CPD__45__12365_c M_3__45__P__45__HYDROXYPYRUVATE_c M_DI__45__H__45__OROTATE_c M_DIHYDRO__45__DIOH__45__BENZOATE_c M_RIBOSE__45__1P_c M_L__45__ALPHA__45__ALANINE_c M_CPD__45__13043_c M_PHOSPHORIBOSYL__45__FORMIMINO__45__AICAR__45__P_c M_PROTON_e M_AMMONIA_c M_INDOLE__45__3__45__GLYCEROL__45__P_c M_P__45__AMINO__45__BENZOATE_c M_CPD__45__8259_c M_GLYCOLALDEHYDE_c M_PHENYL__45__PYRUVATE_c M_HISTIDINOL_c M_NIACINE_c M_N__45__5__45__PHOSPHORIBOSYL__45__ANTHRANILATE_c M_CPD0__45__1108_c M_HIS_c M_3__45__P__45__SERINE_c M_DIAMINO__45__OH__45__PHOSPHORIBOSYLAMINO__45__PYR_c M_GDP__45__D__45__GLUCOSE_c M_OXALO__45__SUCCINATE_c M_NICOTINATE_NUCLEOTIDE_c M_GTP_c M_2__45__KETO__45__3__45__DEOXY__45__6__45__P__45__GLUCONATE_c M_SER_c M_ACET_c M_PPI_c M_GLT_c M_NICOTINAMIDE_RIBOSE_c M_FRUCTOSE__45__16__45__DIPHOSPHATE_c M_ADENINE_c M_CPD__45__62_c M_L__45__ASPARTATE__45__SEMIALDEHYDE_c M_ALPHA__45__D__45__MANNOSYL__45__3__45__PHOSPHOGLYCERATE_c M_TREHALOSE__45__6P_c M_CU__43__2_c M_DAMP_c M_NITRATE_c M_3OH__45__4P__45__OH__45__ALPHA__45__KETOBUTYRATE_c M_XANTHOSINE__45__5__45__PHOSPHATE_c M_CPD0__45__2461_c M_GLN_c M_CPD__45__18118_c M_CARBAMATE_c M_D__45__6__45__P__45__GLUCONO__45__DELTA__45__LACTONE_c M_1__45__L__45__MYO__45__INOSITOL__45__1__45__P_c M_ERYTHRONATE__45__4P_c M_ERYTHROSE__45__4P_c M_4__45__AMINO__45__4__45__DEOXYCHORISMATE_c M_CO3_c M_MYO__45__INOSITOL_c M_D__45__SEDOHEPTULOSE__45__7__45__P_c M_CPD__45__22307_c M_D__45__BETA__45__D__45__HEPTOSE__45__1__45__P_c M_ANTHRANILATE_c M_SULFATE_c M_DGDP_c M_CPD0__45__2472_c M_5__45__P__45__BETA__45__D__45__RIBOSYL__45__AMINE_c M_ENOL__45__PHENYLPYRUVATE_c M_MN__43__2_c M_HISTIDINAL_c M_CPD__45__69_c M_CPD0__45__2101_c M_H2CO3_c M_XTP_c M_SHIKIMATE_c M_CARBAMOYL__45__P_c M_2__45__3__45__DIHYDROXYBENZOATE_c M_NADPH_c M_P3I_c M_L__45__RIBULOSE__45__5__45__P_c M_CPD__45__12377_c M_PHOSPHORIBOSYL__45__ATP_c M_OH_c M_GUANINE_c M_CPD0__45__2474_c M_3__45__DEOXY__45__D__45__ARABINO__45__HEPTULOSONATE__45__7__45__P_c M_FE__43__3_c M_CARBAMYUL__45__L__45__ASPARTATE_c M_NA__43___c M_CPD__45__18238_c M_DIHYDRO__45__NEO__45__PTERIN_c M_CPD__45__16015_c M_SHIKIMATE__45__5P_c M_PHOSPHO__45__ENOL__45__PYRUVATE_c M_TREHALOSE_c M_FRUCTOSE__45__6P_c M_PHOSPHORIBULOSYL__45__FORMIMINO__45__AICAR__45__P_c M_CPD0__45__1699_c M_DPG_c M_L__45__DELTA1__45__PYRROLINE_5__45__CARBOXYLATE_c M_DGMP_c M_INDOLE_ACETALDEHYDE_c M_GLYCEROL__45__3P_c M_D__45__RIBULOSE__45__15__45__P2_c M_5__45__OXOPROLINE_c M_ADP__45__D__45__GLYCERO__45__D__45__MANNO__45__HEPTOSE_c M_OXYGEN__45__MOLECULE_c M_CU__43___e M_D__45__BETA__45__D__45__HEPTOSE__45__17__45__DIPHOSPHATE_c M_7__45__AMINOMETHYL__45__7__45__DEAZAGUANINE_c M_BETA__45__D__45__FRUCTOSE_c M_VAL_c M_D__45__ALANINE_c M_GLC__45__1__45__P_c M_DEAMIDO__45__NAD_c M_CPD__45__3_c M_L__45__DI__45__GMP_c M_D__45__ALA__45__D__45__ALA_c M_CARBOXYPHENYLAMINO__45__DEOXYRIBULOSE__45__P_c M_GUANOSINE__45__5DP__45__3DP_c M_DIHYDRONEOPTERIN__45__P3_c M_ATP_c M_RIBULOSE__45__5P_c M_KDO_c M_GMP_c M_DADP_c M_DGTP_c M_GDP__45__MANNOSE_c M_CPD__45__470_c M_N__45__23__45__DIHYDROXYBENZOYL__45__L__45__SERINE_c M_CPD__45__9923_c M_D__45__glucopyranose__45__6__45__phosphate_c M_ALPHA__45__GLC__45__6__45__P_c M_D__45__GLUCOSAMINE__45__6__45__P_c M_DEOXYINOSINE_c M_GLYCERATE_c M_GLC__45__6__45__P_c M_PROTON_c M_CAMP_c M_FAD_c M_MANNOSE_c M_PAPS_c M_NIACINAMIDE_c M_L__45__LACTATE_c M_CPD__45__302_c M_Glucopyranose_c M_Pi_c M_CIS__45__ACONITATE_c M_CPD__45__2961_c M_FERRIC__45__ENTEROBACTIN__45__COMPLEX_c M_METHYL__45__GLYOXAL_c M_SUC_c M_NMNH_c M_L__45__BETA__45__ASPARTYL__45__P_c M_GDP__45__TP_c M_ZN__43__2_c M_FUM_c M_DIHYDRONEOPTERIN__45__P_c M_GLUCONATE_c M_L__45__GLUTAMATE__45__5__45__P_c M_HYPOXANTHINE_c M_ADENOSINE_c M_D__45__4__45__HYDROXY__45__2__45__KETO__45__GLUTARATE_c M_B__45__ALANINE_c M_3__45__ENOLPYRUVYL__45__SHIKIMATE__45__5P_c M_AICAR_c M_N__45__FORMIMINO__45__L__45__GLUTAMATE_c M_FMNH2_c M_CO__43__2_c M_OH__45__PYR_c M_CPD__45__15979_c M_PREPHENATE_c M_ADENOSYL__45__P4_c M_D__45__LACTATE_c M_CPD__45__407_c M_PHE_c M_2__45__O__45__ALPHA__45__MANNOSYL__45__D__45__GLYCERATE_c M_CPD__45__15382_c M_G3P_c M_ASN_c M_FORMYL__45__ISOGLUTAMINE_c M_CL__45___c M_NICOTINAMIDE_NUCLEOTIDE_c M_2__45__AMINOACRYLATE_c M_GDP_c M_GAP_c M_DIHYDROFOLATE_c M_2__45__ACETO__45__LACTATE_c M_CPD__45__10330_c M_QUINOLINATE_c M_XANTHINE_c M_3__45__DEHYDRO__45__SHIKIMATE_c M_K__43___c M_CPD__45__448_c M_D__45__ALPHABETA__45__D__45__HEPTOSE__45__7__45__PHOSPHATE_c M_CPD__45__1086_c M_DIHYDROXY__45__ACETONE__45__PHOSPHATE_c M_D__45__GLT_c M_OXALACETIC_ACID_c M_DIHYDROPTERIN__45__CH2OH__45__PP_c M_PRPP_c M_7__45__8__45__DIHYDROPTEROATE_c M_2__45__C__45__METHYL__45__D__45__ERYTHRITOL__45__4__45__PHOSPHATE_c M_HYDROGEN__45__PEROXIDE_c M_TRP_c M_DIHYDROMONAPTERIN__45__TRIPHOSPHATE_c M_MALONATE__45__S__45__ALD_c M_CPD__45__15358_c M_CPD__45__10353_c M_HCO3_c M_D__45__Ribofuranose_c M_RIBOFLAVIN_c M_CPD__45__316_c M_PHOSPHORIBOSYL__45__AMP_c M_P__45__HYDROXY__45__PHENYLPYRUVATE_c M_CPD__45__15317_c M_L__45__HISTIDINOL__45__P_c M_CPD__45__12366_c M_AMINO__45__OH__45__HYDROXYMETHYL__45__DIHYDROPTERIDINE_c M_DEOXYXYLULOSE__45__5P_c M_D__45__Ribopyranose_c M_CARBON__45__DIOXIDE_c M_D__45__ERYTHRO__45__IMIDAZOLE__45__GLYCEROL__45__P_c M_L__45__GLUTAMATE_GAMMA__45__SEMIALDEHYDE_c M_3__45__HYDROXY__45__PROPIONATE_c
        intersection of scope 123
        union of scope 325
        max metabolites in scope 321
        min metabolites in scope 127
        average number of metabolites in scope 224.00 (+/- 137.18)
        --- Indiv scopes runtime 0.88 seconds ---
        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/metage2metabo/metage2metabo/workflow_data/out/community_analysis/miscoto_4j9r_2bh.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in out/community_analysis/comm_scopes.json
        --- Community scope runtime 0.73 seconds ---
        Added value of cooperation over individual metabolism: 33 newly reachable metabolites:
        M_DCTP_c M_CPD__45__19306_c M_DEOXYCYTIDINE_c M_CDP_c M_URACIL_c M_UDP__45__D__45__GALACTO__45__14__45__FURANOSE_c M_5__45__HYDROXY__45__CTP_c M_DCDP_c M_CPD__45__12575_c M_2__45__PHOSPHO__45__4__45__CYTIDINE__45__5__45__DIPHOSPHO__45__2__45__C__45__MET_c M_URIDINE_c M_DEOXYADENOSINE_c M_CPD__45__15158_c M_CYTIDINE_c M_CPD__45__16020_c M_CMP_c M_4__45__CYTIDINE__45__5__45__DIPHOSPHO__45__2__45__C_c M_CTP_c M_DEOXYURIDINE_c M_UDP__45__GLUCURONATE_c M_OROTATE_c M_CYTOSINE_c M_UDP_c M_UTP_c M_2C__45__METH__45__D__45__ERYTHRITOL__45__CYCLODIPHOSPHATE_c M_DUTP_c M_CPD__45__14553_c M_THF_c M_CMP__45__KDO_c M_OROTIDINE__45__5__45__PHOSPHATE_c M_UMP_c M_DUMP_c M_DCMP_c
        Setting these 33 as targets
        Running minimal community selection
        Community scopes for all metabolic networks available in out/community_analysis/comm_scopes.json
        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities = 2
        GCA_003433675
        GCA_003433665
        ######### Keystone species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Keystone species = 2
        GCA_003433675
        GCA_003433665
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal community enabling the producibility of the target metabolites given as inputs
        Essential symbionts = 2
        GCA_003433675
        GCA_003433665
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal community enabling the producibility of the target metabolites given as inputs
        Alternative symbionts = 0
        --- Mincom runtime 1.02 seconds ---
        --- Total runtime 538.29 seconds ---

* files outputs
    * Numerous files are created in the output_directory, including the logs at the root of the results directory.
    
    .. code ::

        output_directory/
        ├── m2m_workflow.log
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── mincom.json
        │   ├── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        ├── padmet
        │   ├── GCA_003433665.padmet
        │   └── GCA_003433675.padmet
        ├── pgdb
        │   ├── GCA_003433665
        │   │   ├── classes.dat
        │   │   ├── compound-links.dat
        │   │   ├── compounds.dat
        │   │   ├── dnabindsites.dat
        │   │   ├── enzrxns.dat
        │   │   ├── gene-links.dat
        │   │   ├── genes.dat
        │   │   ├── pathway-links.dat
        │   │   ├── pathways.dat
        │   │   ├── promoters.dat
        │   │   ├── protein-features.dat
        │   │   ├── protein-links.dat
        │   │   ├── proteins.dat
        │   │   ├── protligandcplxes.dat
        │   │   ├── pubs.dat
        │   │   ├── reaction-links.dat
        │   │   ├── reactions.dat
        │   │   ├── regulation.dat
        │   │   ├── regulons.dat
        │   │   ├── rnas.dat
        │   │   ├── species.dat
        │   │   ├── terminators.dat
        │   │   └── transunits.dat
        │   └── GCA_003433675
        │       ├── classes.dat
        │       ├── compound-links.dat
        │       ├── compounds.dat
        │       ├── dnabindsites.dat
        │       ├── enzrxns.dat
        │       ├── gene-links.dat
        │       ├── genes.dat
        │       ├── pathway-links.dat
        │       ├── pathways.dat
        │       ├── promoters.dat
        │       ├── protein-features.dat
        │       ├── protein-links.dat
        │       ├── proteins.dat
        │       ├── protligandcplxes.dat
        │       ├── pubs.dat
        │       ├── reaction-links.dat
        │       ├── reactions.dat
        │       ├── regulation.dat
        │       ├── regulons.dat
        │       ├── rnas.dat
        │       ├── species.dat
        │       ├── terminators.dat
        │       └── transunits.dat
        └── recon_stats.tsv
        └── sbml
            ├── GCA_003433665.sbml
            └── GCA_003433675.sbml

    These files are the same as the ones presented in the previous commands: metabolic networks reconstructions (Pathway Tools data, SBML), individual and collective scopes, minimal community selection. 


Including a host in the picture
-------------------------------

It is possible to consider a host in addition to the microbiota for the ``workflow``, ``cscope`` and ``mincom`` commands. **What does it change?**

First note that adding the host in the SBML repository will enable you to get the individual scope for the host. Another solution is to directly use ``menescope`` from the `MeneTools
<https://github.com/cfrioux/MeneTools>`_ `Python package <https://pypi.org/project/MeneTools/>`__ on which m2m relies, and that can be used as a standalone tool.

Then back to the effect of the host in the other commands.

* For ``cscope`` and ``addedvalue``, the host metabolism will be taken into account. That is to say that it will be considered as a member of the community. Among the newly producible targets, some will be exclusive to the host metabolism. This is not displayed in the standard output of the software but can be retrieved in the json file output under the `"comhost_scope"` key of the dictionary. 

* For ``mincom``, the host will always be considered in the community. This means that the selected bacteria need to be associated to the host in order to ensure the producibility of all the targets. Therefore, if the minimal community computed for 10 targets is of 3 bacteria and that a host was provided, it means that the host + these three bacteria can produce the 10 targets. 

More generally, for more information and analysis on the usage of hosts in addition to the microbiota, we refer the interested user to the `Miscoto
<https://github.com/cfrioux/miscoto>`_ `Python package <https://pypi.org/project/Miscoto/>`__, on which m2m relies. Miscoto can be used as a standalone package for such analyses, with additional options, such as the identification of putative exchanges among the minimal communities. 
