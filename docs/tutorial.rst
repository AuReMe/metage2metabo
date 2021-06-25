============
m2m Tutorial
============
Test data is avaible in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/test>`__.
It contains enough data to run the different subcommands.

Outputs and logs
----------------

By default, ``m2m`` writes into the console (stdout) and into a log file located at the root of the results directory and named after the subcommand that was executed. The option ``-q`` can be given to any ``m2m`` subcommand to write in the log file and to stdout only the warnings, errors and critical issues, together with the path to the log file.

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
--pwt-xml        extract xml from Pathway Tools instead of using padmet to create sbml

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

In addition, the genbank files (``.gbk`` or ``.gbff``) also follow the strict requirements of Pathway Tools.
Please go check the `documentation of mpwt <https://github.com/AuReMe/mpwt#genbank>`__ for more details, especially if you get errors from Pathway Tools at the ``recon`` step.

.. code:: sh

    m2m recon -g workflow_genomes -o output_directory -c cpu_number [--clean] [--orphan] [-p] [-l sbml_level] [--pwt-xml]

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
    * In ``output_directory/pgdb``, the .dat files of Pathway Tools. The corresponding SBMLs are in ``output_directory/sbml``. The structure of the output directory after this ``recon`` command is shown below :

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

By using the ``--pwt-xml``, m2m will extract the xml files created by MetaFlux from the PGDBs created by Pathway Tools. This will modify the files stored in the ``pgdb`` folder:

::

    output_directory/
    ├── m2m_recon.log
    ├── pgdb
    │   ├── GCA_003433665.xml
    │   └── GCA_003433675.xml
    └── recon_stats.tsv
    └── sbml
        ├── GCA_003433665.sbml
        └── GCA_003433675.sbml

m2m iscope, cscope and addedvalue
---------------------------------
The three subcommands require metabolic networks in the SBML format. Some metabolic networks are available as a compressed archive in `metabolic_data`. Uncompress the file and the directory can be fed to the subcommands. These commands also require a seeds file comprising the metabolic compounds available to assess reachability/producibility in the models. This seeds file needs to be in SBML format. You can use the one in the ``metabolic_data`` directory.

iscope
*******

It uses the following mandatory inputs (run ``m2m iscope --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional argument

-q                     quiet mode
-c int           number of CPU for multi-processing

.. code:: sh

    m2m iscope -n toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory/

* standard output
    .. code:: 

        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory//indiv_scopes/indiv_scopes.json
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
        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory//indiv_scopes/rev_iscope.json and as a matrix in output_directory//indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 9.07 seconds ---

        --- Total runtime 9.08 seconds ---
        --- Logs written in output_directory//m2m_iscope.log ---

    These results mean that 135 metabolites can be reached by all organisms. When gathering reachable metabolites for all organisms, the union consists of 625 metabolites (including the seeds). Finally metrics show the min, max and arithmetic mean number of compounds in all scopes.
* files outputs
    * In ``output_directory/indiv_scopes/indiv_scopes.json``: a json file that can be easily loaded as a dictionary (or humanly read as it it) that contains the set of reachable metabolites for each organism. /!\\ Warning: the seeds are included in the scopes, hence they will never be empty. Logs are written in ``output_directory/m2m_iscope.log``.

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
-t file          Optional targets for metabolic analysis, if not used
                 metage2metabo will use the addedvalue of the community
-q                     quiet mode

.. code:: sh

    m2m cscope -n toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory/

* standard output
    .. code::

        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/programs/metage2metabo/test/output_directory/community_analysis/miscoto_96owqje2.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory//community_analysis/comm_scopes.json
        --- Community scope runtime 4.33 seconds ---


        651 metabolites (excluding the seeds) reachable by the whole community/microbiota: 

        ...
        --- Total runtime 4.34 seconds ---
        --- Logs written in output_directory//m2m_cscope.log ---

    651 metabolites are reachable by the microbiota. This does not include the seeds. The list of metabolites is given in output. 
* files outputs
    * In addition to the logs at the root of the results directory, a json file with the results is created in ``output_directory/community_analysis/indiv_scopes.json``.
    * To screen the putative redundancy of metabolite producibility predicted for the genomes, we also provide `rev_iscope.tsv` and `rev_iscope.json` that reverse the information from `indiv_scopes.json`. This means that if org1 produces A and B, org2 produces B and C, `indiv_scopes.json` will describe the following: {'org1': ['A', 'B'], 'org2: ['B', 'C']}. `reverse_scope.json` will contain {'A': ['org1'], 'B': ['org1', 'org2'], 'C': ['org2']}, and `reverse_scope.tsv` will contain the same information as a matrix. 

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
        Individual scopes for all metabolic networks available in output_directory//indiv_scopes/indiv_scopes.json
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
        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory//indiv_scopes/rev_iscope.json and as a matrix in output_directory//indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 8.43 seconds ---

        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/programs/metage2metabo/test/output_directory/community_analysis/miscoto_rwvbc87b.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory//community_analysis/comm_scopes.json
        --- Community scope runtime 4.33 seconds ---


        651 metabolites (excluding the seeds) reachable by the whole community/microbiota: 

        ...

        Added value of cooperation over individual metabolism: 119 newly reachable metabolites: 

        ...

        Added-value of cooperation written in output_directory//community_analysis/addedvalue.json
        Target file created with the addedvalue targets in: output_directory//community_analysis/targets.sbml
        --- Total runtime 12.78 seconds ---
        --- Logs written in output_directory//m2m_addedvalue.log ---

    As you can see, the individual and community scopes are run again. In addition to the previous outputs, the union of all individual scopes and the community scopes are printed. Finally, the difference between the two sets, that is to say the metabolites that can only be produced collectively (i.e. by at least two bacteria cooperating) is displayed. Here it consists of 119 metabolites. 
* files outputs
    * A targets SBML file is generated. It can be used with ``m2m mincom`` . Newly producible metabolites are written in a json file. The json files associated to ``iscope`` and ``cscope`` are also produced.

    ::

        output_directory/
        ├── m2m_addedvalue.log
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── addedvalue.json
        │   └── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv

Optional: create the seeds SBML file
*************************************
To create a seeds file starting from a list of metabolic identifiers (matching identifiers of compounds of the organisms metabolic networks), you can use the ``m2m seeds`` command:

.. code:: sh

    m2m seeds --metabolites metabolites_file.txt -o output/directory

The resulting seeds file will be created in ``output/directory/seeds.sbml``.

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

m2m mincom
----------
``m2m mincom`` requires an additional target file that is available in `metabolic_data` or can be generated by ``m2m addedvalue`` in which case it will be stored in ``result_directory/community_analysis/targets.sbml``

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
        Created instance in /shared/programs/metage2metabo/test/output_directory/community_analysis/miscoto_fsboc3q7.lp
        Running minimal community selection

        In the initial and minimal communities 120 targets are producible and 0 remain unproducible.

        120 producible targets:
        ...

        0 still unproducible targets:


        Minimal communities are available in output_directory//community_analysis/mincom.json 

        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities => 13

        GCA_003437885
        GCA_003437715
        GCA_003437905
        GCA_003437665
        GCA_003437255
        GCA_003437295
        GCA_003437815
        GCA_003437945
        GCA_003438055
        GCA_003437195
        GCA_003437375
        GCA_003437055
        GCA_003437595
        ######### Key species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Number of key species => 17

        GCA_003437885
        GCA_003437715
        GCA_003437905
        GCA_003437665
        GCA_003437295
        GCA_003438055
        GCA_003437195
        GCA_003437175
        GCA_003437345
        GCA_003437055
        GCA_003437595
        GCA_003437785
        GCA_003437255
        GCA_003437815
        GCA_003437945
        GCA_003437375
        GCA_003437325
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs
        Number of essential symbionts => 12

        GCA_003437885
        GCA_003437715
        GCA_003437905
        GCA_003437665
        GCA_003437255
        GCA_003437295
        GCA_003437815
        GCA_003438055
        GCA_003437195
        GCA_003437375
        GCA_003437055
        GCA_003437595
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs
        Number of alternative symbionts => 5

        GCA_003437785
        GCA_003437945
        GCA_003437345
        GCA_003437175
        GCA_003437325

        --- Mincom runtime 3.70 seconds ---

        --- Total runtime 6.50 seconds ---
        --- Logs written in output_directory//m2m_mincom.log ---


    This output gives the result of minimal community selection. It means that for producing the 120 metabolic targets, a minimum of 13 bacteria out of the 17 is required. One example of such minimal community is given. In addition, the whole space of solution is studied. All bacteria (17) occur in at least one minimal community (key species). Finally, the intersection gives the following information: a set of 12 bacteria occurs in each minimal communtity. This means that these 12 bacteria are needed in any case (essential symbionts), and that any of the remaining 5 bacteria (alternative symbionts) can complete the missing function(s).
* files outputs
    * As for other commands, a json file with the results is produced in ``output_directory/community_analysis/comm_scopes.json``, together with logs at the root of the results directory.

m2m metacom
------------
``m2m metacom`` runs all analyses: individual scopes, community scopes, and minimal community selection based on the metabolic added-value of the microbiota.

It uses the following mandatory inputs (run ``m2m metacom --help`` for optional arguments):

-n directory           directory of metabolic networks,
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional arguments:

-m file                host metabolic network SBML file
-t file                Optional targets for metabolic analysis, if not used
                       metage2metabo will use the addedvalue of the community
-q                     quiet mode
-c int           number of CPU for multi-processing

.. code:: sh

    m2m metacom -n metabolic_data/toy_bact -s metabolic_data/seeds_toy.sbml  -o output_directory

* standard output
    .. code::

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
        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory/indiv_scopes/rev_iscope.json and as a matrix in output_directory/indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 9.77 seconds ---

        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/programs/metage2metabo/test/output_directory/community_analysis/miscoto_wkdkeazl.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        --- Community scope runtime 5.84 seconds ---


        Added value of cooperation over individual metabolism: 119 newly reachable metabolites: 

        ...


        Added-value of cooperation written in output_directory/community_analysis/addedvalue.json

        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml
        Setting 119 compounds as targets 

        Running minimal community selection

        In the initial and minimal communities 119 targets are producible and 0 remain unproducible.

        119 producible targets:
        ...

        0 still unproducible targets:


        Minimal communities are available in output_directory/community_analysis/mincom.json 

        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities => 13

        GCA_003437255
        GCA_003437885
        GCA_003437815
        GCA_003437375
        GCA_003437295
        GCA_003437715
        GCA_003437665
        GCA_003438055
        GCA_003437195
        GCA_003437905
        GCA_003437595
        GCA_003437055
        GCA_003437945
        ######### Key species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Number of key species => 17

        GCA_003437785
        GCA_003437885
        GCA_003437055
        GCA_003437345
        GCA_003437665
        GCA_003437195
        GCA_003437905
        GCA_003437175
        GCA_003437595
        GCA_003437325
        GCA_003437815
        GCA_003437375
        GCA_003437295
        GCA_003437715
        GCA_003437255
        GCA_003438055
        GCA_003437945
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs
        Number of essential symbionts => 12

        GCA_003437255
        GCA_003437885
        GCA_003437815
        GCA_003437295
        GCA_003437375
        GCA_003437715
        GCA_003437665
        GCA_003438055
        GCA_003437195
        GCA_003437905
        GCA_003437595
        GCA_003437055
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs
        Number of alternative symbionts => 5

        GCA_003437345
        GCA_003437945
        GCA_003437175
        GCA_003437325
        GCA_003437785

        --- Mincom runtime 4.34 seconds ---

        Targets producibility are available at output_directory/producibility_targets.json
        --- Total runtime 20.01 seconds ---
        --- Logs written in output_directory/m2m_metacom.log ---

* files outputs
    * Files are created in the output_directory: the logs, json files with the results, targets in SBML.

    .. code ::

        output_directory/
        ├── m2m_metacom.log
        ├── producibility_targets.json
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── mincom.json
        │   ├── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv

m2m workflow and m2m test
-------------------------
``m2m workflow`` starts from metabolic network reconstruction and runs all analyses: individual scopes, community scopes, and minimal community selection based on the metabolic added-value of the microbiota.

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
-t file          Optional targets for metabolic analysis, if not used
                 metage2metabo will use the addedvalue of the community
-q               quiet mode
--pwt-xml        extract xml from Pathway Tools instead of using padmet to create sbml

You can run the workflow analysis with the two genbanks files available in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/metage2metabo>`__ (`workflow_data`). Two genomes are available in the compressed archive `workflow_genomes.tar.gz`. The archive has to be uncompressed before testing.

.. code:: sh

    m2m workflow -g workflow_genomes -s workflow_data/seeds_workflow.sbml -o output_directory/

Or you can run the test argument (which use the same data): ``m2m test``.

Which uses the following mandatory inputs (run ``m2m test --help`` for optional arguments):

-o directory           output directory path

Optional arguments:

-q               quiet mode
-c int           cpu number for multi-processing

.. code:: sh

    m2m test -o output_directory

* standard outputs

    .. code ::

        Uncompressing test data to output_directory
        Launching workflow on test data
        ######### Running metabolic network reconstruction with Pathway Tools #########
        ~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~
        Checking inputs for GCA_003433675: no missing files.
        Checking inputs for GCA_003433665: no missing files.
        ----------End of creation of input data from Genbank/GFF/PF: 0.02s----------
        ~~~~~~~~~~Inference on the data~~~~~~~~~~
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho output_directory/workflow_genomes/GCA_003433665/
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho output_directory/workflow_genomes/GCA_003433675/
        ~~~~~~~~~~Check inference~~~~~~~~~~
        No log directory, it will be created.

        2 builds have passed!

        ----------End of PathoLogic inference: 403.90s----------
        ~~~~~~~~~~Creation of the .dat files~~~~~~~~~~
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load output_directory/workflow_genomes/GCA_003433675/dat_creation.lisp
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load output_directory/workflow_genomes/GCA_003433665/dat_creation.lisp
        ~~~~~~~~~~Check .dat~~~~~~~~~~
        gca_003433675cyc: 23 out of 23 dat files created.
        gca_003433665cyc: 23 out of 23 dat files created.
        ----------End of dat files creation: 163.21s----------
        ~~~~~~~~~~End of Pathway Tools~~~~~~~~~~
        ~~~~~~~~~~Moving result files~~~~~~~~~~
        ----------End of moving fimes: 0.12s----------
        ----------mpwt has finished in 567.29s! Thank you for using it.
        ######### Creating SBML files #########
        ######### Stats GSMN reconstruction #########
        Number of genomes: 2
        Number of reactions in all GSMN: 2026
        Number of compounds in all GSMN: 2095
        Average reactions per GSMN: 1437.00(+/- 678.82)
        Average compounds per GSMN: 1560.00(+/- 615.18)
        Average genes per GSMN: 893.00(+/- 475.18)
        Average pathways per GSMN: 257.00(+/- 134.35)
        Percentage of reactions associated with genes: 79.90(+/- 3.20)
        --- Recon runtime 574.26 seconds ---

        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json
        2 metabolic models considered.

        123 metabolites in core reachable by all organisms (intersection) 

        ...

        325 metabolites reachable by individual organisms altogether (union), among which 26 seeds (growth medium) 

        ...

        intersection of scope 123
        union of scope 325
        max metabolites in scope 321
        min metabolites in scope 127
        average number of metabolites in scope 224.00 (+/- 137.18)
        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory/indiv_scopes/rev_iscope.json and as a matrix in output_directory/indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 1.21 seconds ---

        ######### Creating metabolic instance for the whole community #########
        Created instance in /shared/programs/metage2metabo/test/output_directory/community_analysis/miscoto_17f5ygw7.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        --- Community scope runtime 0.85 seconds ---


        Added value of cooperation over individual metabolism: 33 newly reachable metabolites: 

        ...


        Added-value of cooperation written in output_directory/community_analysis/addedvalue.json

        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml
        Setting 33 compounds as targets 

        Running minimal community selection

        In the initial and minimal communities 33 targets are producible and 0 remain unproducible.

        33 producible targets:
        ...

        0 still unproducible targets:


        Minimal communities are available in output_directory/community_analysis/mincom.json 

        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities => 2

        GCA_003433665
        GCA_003433675
        ######### Key species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Number of key species => 2

        GCA_003433665
        GCA_003433675
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs
        Number of essential symbionts => 2

        GCA_003433665
        GCA_003433675
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs
        Number of alternative symbionts => 0



        --- Mincom runtime 1.36 seconds ---

        Targets producibility are available at output_directory/producibility_targets.json
        --- Total runtime 577.84 seconds ---
        --- Logs written in output_directory/m2m_test.log ---


* files outputs
    * Numerous files are created in the `output_directory`, including the logs at the root of the results directory.
    
    .. code ::

        output_directory/
        ├── m2m_workflow.log
        ├── producibility_targets.json
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── mincom.json
        │   ├── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
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
