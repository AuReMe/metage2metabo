============
m2m Tutorial
============
Test data is avaible in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/main/test>`__.
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

This structure can be observed in the `workflow_genomes.tar.gz` file in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/main/metage2metabo/workflow_data>`__.

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

        ###############################################
        #                                             #
        #       Metabolic network reconstruction      #
        #                                             #
        ###############################################

        ######### Running metabolic network reconstruction with Pathway Tools #########
        ---------- Launching mpwt ----------
        |Input Check|GCA_003433665| No missing files
        |PathoLogic|GCA_003433665| pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_003433665 -dump-flat-files-biopax
        |Input Check|GCA_003433675| No missing files
        |PathoLogic|GCA_003433675| pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_003433675 -dump-flat-files-biopax
        |Output Check|workflow_genomes/GCA_003433675| 23 out of 23 dat files created.
        |Moving output files|GCA_003433675| 
        |Output Check|workflow_genomes/GCA_003433665| 23 out of 23 dat files created.
        |Moving output files|GCA_003433665| 
        |Output Check| 2 on 2 builds have passed!
        -------------- Checking mpwt runs --------------
        All runs are successful.
        -------------- mpwt has finished in 251.79s! Thank you for using it. --------------
        ######### Creating SBML files #########
        ######### Stats GSMN reconstruction #########
        Number of genomes: 2
        Number of reactions in all GSMN: 2432
        Number of compounds in all GSMN: 2387
        Average reactions per GSMN: 1743.00(+/- 780.65)
        Average compounds per GSMN: 1776.00(+/- 695.79)
        Average genes per GSMN: 932.00(+/- 504.87)
        Percentage of reactions associated with genes: 67.95(+/- 4.90)
        --- Recon runtime 256.75 seconds ---

        PGDB created in output_directory/pgdb
        SBML files created in output_directory/sbml
        --- Total runtime 256.75 seconds ---
        --- Logs written in output_directory/m2m_recon.log ---


        The output shows that PGDB are created with Pathway Tools. Then the .dat files are extracted and used to build SBML files of the metabolic models.
* files outputs
    * In ``output_directory/pgdb``, the .dat files of Pathway Tools. The corresponding SBMLs are in ``output_directory/sbml``. The structure of the output directory after this ``recon`` command is shown below :

    ::

        output_directory/
        ├── m2m_metadata.json
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
    ├── m2m_metadata.json
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
******

It uses the following mandatory inputs (run ``m2m iscope --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

Optional argument

-q                     quiet mode
-c int           number of CPU for multi-processing

.. code:: sh

    m2m iscope -n toy_bact -s seeds_toy.sbml -o output_directory/

* standard output

    .. code:: 

        ###############################################
        #                                             #
        #       Individual metabolic potentials       #
        #                                             #
        ###############################################


        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json. The scopes have been filtered a way that if a seed is in a scope, it means the corresponding species is predicted to be able to produce it.

        Information regarding the producibility of seeds, and the possible absence of seeds in some metabolic networks is stored in output_directory/indiv_scopes/seeds_in_indiv_scopes.json.

        17 metabolic models considered.

        50 metabolites in core reachable by all organisms (intersection)
        
        ...

        576 metabolites reachable by individual organisms altogether (union), among which 44 metabolites that are also part of the seeds (growth medium)

        ...

        Summary:
        - intersection of scopes 50
        - union of scopes 576
        - max metabolites in scopes 422
        - min metabolites in scopes 116
        - average number of metabolites in scopes 239.06 (+/- 89.51)

        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory/indiv_scopes/rev_iscope.json and as a matrix in output_directory/indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 21.46 seconds ---

        --- Total runtime 21.47 seconds ---
        --- Logs written in output_directory/m2m_iscope.log --
    

    These results mean that 50 metabolites can be reached by all organisms. When gathering reachable metabolites for all organisms, the union consists of 576 metabolites. Some of the reachable metabolites can also be part of the seeds, meaning that there would be a possibility to renew the reservoir of seed molecules by some species. Finally metrics show the min, max and arithmetic mean number of compounds in all scopes.

* files outputs

    * In ``output_directory/indiv_scopes/indiv_scopes.json``: a json file that can be easily loaded as a dictionary (or humanly read as it it) that contains the set of reachable metabolites for each organism.
    * A file expliciting the producibility of seeds ``output_directory/indiv_scopes/seeds_in_indiv_scopes.json`` is also available: it additionally lists the seeds that are absent from the networks.
    * Two more files present the scopes from the focus of metabolites ``output_directory/indiv_scopes/rev_iscopes.json`` and a matrix summarising the producibility of molecules by species ``output_directory/indiv_scopes/rev_iscopes.tsv``. `rev_iscope.tsv` and `rev_iscope.json` that reverse the information from `indiv_scopes.json`. This means that if org1 produces A and B, org2 produces B and C, `indiv_scopes.json` will describe the following: {'org1': ['A', 'B'], 'org2: ['B', 'C']}. `reverse_scope.json` will contain {'A': ['org1'], 'B': ['org1', 'org2'], 'C': ['org2']}, and `reverse_scope.tsv` will contain the same information as a matrix.
    * Logs are written in ``output_directory/m2m_iscope.log`` and metadata containing package versions in ``output_directory/m2m_metadata.json``.

cscope
******

It uses the following mandatory inputs (run ``m2m cscope --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file                targets SBML file
-o directory           output directory for results

Optional arguments:

-m file                host metabolic network SBML file
-t file          Optional targets for metabolic analysis, if not used
                 metage2metabo will use the addedvalue of the community
-q                     quiet mode

.. code:: sh

    m2m cscope -n toy_bact -s seeds_toy.sbml -o output_directory/

* standard output

    .. code::

        ###############################################
        #                                             #
        #    Metabolic potential of the community     #
        #                                             #
        ###############################################

        ######### Creating metabolic instance for the whole community #########
        Created temporary instance file in ../metage2metabo/test/metabolic_data/output_directory/community_analysis/miscoto_9ihtb055.lp
        Running whole-community metabolic scopes...
        Community scope for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        Contributions of microbes to community scope available in output_directory/community_analysis/contributions_of_microbes.json.


        Number of metabolites producible in community: 698.

        Reverse community scopes for all metabolic networks available in output_directory/community_analysis/rev_cscope.json and output_directory/community_analysis/rev_cscope.tsv. They higlight the producibility of metabolites by species in the community.

        --- Community scope runtime 5.41 seconds ---
        ...
        --- Total runtime 5.42 seconds ---
        --- Logs written in output_directory/m2m_cscope.log ---


698 metabolites (excluding the seeds) reachable by the whole community/microbiota:

* files outputs
    * In addition to the logs at the root of the results directory, a json file with the results is created in ``output_directory/community_analysis/comm_scopes.json``. It lists all molecules reachable by the community, taking into account the interactions occurring among community members.
    * A file details the roles of community members in the production of metabolites: which microbes possess the reactions that produce the metabolites. This file is ``output_directory/community_analysis/contributions_of_microbes.json``. It also recapitulates the compounds producible by species individually versus in community, and highlights the newly producible compounds in community, per symbiont. 
    * As for the individual scopes, the redundancy of metabolite producibility is described in ``output_directory/community_analysis/rev_cscope.json`` and ``output_directory/community_analysis/rev_cscope.tsv``.
    * Logs are written in ``output_directory/m2m_cscope.log`` and metadata containing package versions in ``output_directory/m2m_metadata.json``.

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

    m2m addedvalue -n toy_bact -s seeds_toy.sbml -o output_directory/

* standard output
    .. code::

        ###############################################
        #                                             #
        #       Individual metabolic potentials       #
        #                                             #
        ###############################################


        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json. The scopes have been filtered a way that if a seed is in a scope, it means the corresponding species is predicted to be able to produce it.

        Information regarding the producibility of seeds, and the possible absence of seeds in some metabolic networks is stored in output_directory/indiv_scopes/seeds_in_indiv_scopes.json.

        17 metabolic models considered.

        50 metabolites in core reachable by all organisms (intersection)
        
        ...

        576 metabolites reachable by individual organisms altogether (union), among which 44 metabolites that are also part of the seeds (growth medium)

        ...

        Summary:
        - intersection of scopes 50
        - union of scopes 576
        - max metabolites in scopes 422
        - min metabolites in scopes 116
        - average number of metabolites in scopes 239.06 (+/- 89.51)

        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory/indiv_scopes/rev_iscope.json and as a matrix in output_directory/indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 21.46 seconds ---

        ###############################################
        #                                             #
        #    Metabolic potential of the community     #
        #                                             #
        ###############################################

        ######### Creating metabolic instance for the whole community #########
        Created temporary instance file in ../metage2metabo/test/metabolic_data/output_directory/community_analysis/miscoto_9ihtb055.lp
        Running whole-community metabolic scopes...
        Community scope for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        Contributions of microbes to community scope available in output_directory/community_analysis/contributions_of_microbes.json.


        Number of metabolites producible in community: 698.

        Reverse community scopes for all metabolic networks available in output_directory/community_analysis/rev_cscope.json and output_directory/community_analysis/rev_cscope.tsv. They higlight the producibility of metabolites by species in the community.

        --- Community scope runtime 5.41 seconds ---
        ...

        ###############################################
        #                                             #
        #    Added-value of metabolic interactions    #
        #                                             #
        ###############################################


        Added value of cooperation over individual metabolism: 122 newly reachable metabolites:
        ...
        Added-value of cooperation written in output_directory/community_analysis/addedvalue.json
        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml
        --- Total runtime 27.74 seconds ---
        --- Logs written in output_directory/m2m_addedvalue.log ---

    As you can see, the individual and community scopes are run again. In addition to the previous outputs, the union of all individual scopes and the community scopes are printed. Finally, the difference between the two sets, that is to say the metabolites that can only be produced collectively (i.e. by at least two bacteria cooperating) is displayed. Here it consists of 119 metabolites. 
* files outputs
    * A targets SBML file is generated. It can be used with ``m2m mincom`` . Newly producible metabolites are written in a json file. The json files associated to ``iscope`` and ``cscope`` are also produced.

    ::

        output_directory/
        ├── m2m_addedvalue.log
        ├── m2m_metadata.json
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── addedvalue.json
        │   ├── contributions_of_microbes.json
        │   ├── rev_cscope.json
        │   ├── rev_cscope.tsv
        │   └── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
        │   └── seeds_in_indiv_scopes.json

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

    m2m mincom -n toy_bact -s seeds_toy.sbml -t targets_toy.sbml -o output_directory/

* standard output
    .. code::

        ###############################################
        #                                             #
        #         Minimal community selection         #
        #                                             #
        ###############################################

        WARNING: The following seeds are among the targets: {'M_MANNITOL_c'}. They will not be considered as targets during the computation of minimal communities: they will be considered as already reachable according to the network expansion definition.

        Running minimal community selection
        /Users/cfrioux/.pyenv/versions/metage2metabo/lib/python3.10/site-packages/miscoto/encodings/community_soup.lp

        In the initial and minimal communities 120 targets are producible and 0 remain unproducible.

        120 producible targets:
        ...

        0 still unproducible targets:


        Minimal communities are available in output_directory/community_analysis/mincom.json

        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities => 13

        GCA_003437055
        GCA_003437595
        GCA_003437295
        GCA_003437345
        GCA_003437715
        GCA_003437815
        GCA_003437905
        GCA_003437375
        GCA_003437195
        GCA_003438055
        GCA_003437885
        GCA_003437665
        GCA_003437255
        ######### Key species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Number of key species => 17

        GCA_003437055
        GCA_003437325
        GCA_003437595
        GCA_003437345
        GCA_003437715
        GCA_003437905
        GCA_003437945
        GCA_003438055
        GCA_003437255
        GCA_003437295
        GCA_003437785
        GCA_003437815
        GCA_003437175
        GCA_003437375
        GCA_003437195
        GCA_003437885
        GCA_003437665
        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs
        Number of essential symbionts => 12

        GCA_003437055
        GCA_003437595
        GCA_003437295
        GCA_003437715
        GCA_003437815
        GCA_003437905
        GCA_003437375
        GCA_003437195
        GCA_003438055
        GCA_003437885
        GCA_003437665
        GCA_003437255
        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs
        Number of alternative symbionts => 5

        GCA_003437345
        GCA_003437325
        GCA_003437945
        GCA_003437785
        GCA_003437175

        --- Mincom runtime 5.61 seconds ---

        --- Total runtime 7.72 seconds ---
        --- Logs written in output_directory/m2m_mincom.log ---

    This output gives the result of minimal community selection. It means that for producing the 120 metabolic targets, a minimum of 13 bacteria out of the 17 is required. One example of such minimal community is given. In addition, the whole space of solution is studied. All bacteria (17) occur in at least one minimal community (key species). Finally, the intersection gives the following information: a set of 12 bacteria occurs in each minimal communtity. This means that these 12 bacteria are needed in any case (essential symbionts), and that any of the remaining 5 bacteria (alternative symbionts) can complete the missing function(s).
* files outputs
    * As for other commands, a json file with the results is produced in ``output_directory/community_analysis/mincom.json``, together with logs at the root of the results directory.

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
--target-com-scope           Instead of the addedvalue, use the community scope as targets for mincom

.. code:: sh

    m2m metacom -n metabolic_data/toy_bact -s metabolic_data/seeds_toy.sbml -o output_directory

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
        ├── m2m_metadata.json
        ├── m2m_mincom.log
        ├── producibility_targets.json
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── contributions_of_microbes.json
        │   ├── mincom.json
        │   ├── rev_cscope.json
        │   ├── rev_cscope.tsv
        │   └── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
        │   └── seeds_in_indiv_scopes.json


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
-m file                host metabolic network SBML file
--pwt-xml        extract xml from Pathway Tools instead of using padmet to create sbml
--target-com-scope           Instead of the addedvalue, use the community scope as targets for mincom

You can run the workflow analysis with the two genbanks files available in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/main/metage2metabo>`__ (`workflow_data`). Two genomes are available in the compressed archive `workflow_genomes.tar.gz`. The archive has to be uncompressed before testing.

.. code:: sh

    m2m workflow -g workflow_genomes -s seeds_workflow.sbml -o output_directory/

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

        ###############################################
        #                                             #
        #       Metabolic network reconstruction      #
        #                                             #
        ###############################################

        ######### Running metabolic network reconstruction with Pathway Tools #########
        ---------- Launching mpwt ----------
        |Input Check|GCA_003433665| No missing files
        |PathoLogic|GCA_003433665| pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_003433665 -dump-flat-files-biopax
        |Output Check|workflow_genomes/GCA_003433665| 23 out of 23 dat files created.
        |Moving output files|GCA_003433665| 
        |Input Check|GCA_003433675| No missing files
        |PathoLogic|GCA_003433675| pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_003433675 -dump-flat-files-biopax
        |Output Check|workflow_genomes/GCA_003433675| 23 out of 23 dat files created.
        |Moving output files|GCA_003433675| 
        |Output Check| 2 on 2 builds have passed!
        -------------- Checking mpwt runs --------------
        All runs are successful.
        -------------- mpwt has finished in 415.62s! Thank you for using it. --------------
        ######### Creating SBML files #########
        ######### Stats GSMN reconstruction #########
        Number of genomes: 2
        Number of reactions in all GSMN: 2433
        Number of compounds in all GSMN: 2389
        Average reactions per GSMN: 1743.50(+/- 781.35)
        Average compounds per GSMN: 1777.00(+/- 697.21)
        Average genes per GSMN: 932.00(+/- 504.87)
        Percentage of reactions associated with genes: 67.96(+/- 4.91)
        --- Recon runtime 421.74 seconds ---


        ###############################################
        #                                             #
        #       Individual metabolic potentials       #
        #                                             #
        ###############################################


        Individual scopes for all metabolic networks available in output_directory/indiv_scopes/indiv_scopes.json. The scopes have been filtered a way that if a seed is in a scope, it means the corresponding species is predicted to be able to produce it.

        Information regarding the producibility of seeds, and the possible absence of seeds in some metabolic networks is stored in output_directory/indiv_scopes/seeds_in_indiv_scopes.json.

        2 metabolic models considered.

        143 metabolites in core reachable by all organisms (intersection) 

        ...

        334 metabolites reachable by individual organisms altogether (union), among which 12 metabolites that are also part of the seeds (growth medium) 

        ...

        Summary:
        - intersection of scopes 143
        - union of scopes 334
        - max metabolites in scopes 333
        - min metabolites in scopes 144
        - average number of metabolites in scopes 238.50 (+/- 133.64)

        Analysis of functional redundancy (producers of all metabolites) is computed as a dictionary in output_directory/indiv_scopes/rev_iscope.json and as a matrix in output_directory/indiv_scopes/rev_iscope.tsv.
        --- Indiv scopes runtime 0.83 seconds ---


        ###############################################
        #                                             #
        #    Metabolic potential of the community     #
        #                                             #
        ###############################################

        ######### Creating metabolic instance for the whole community #########
        Created temporary instance file in /shared/Softwares/git/metage2metabo/metage2metabo/workflow_data/output_directory/community_analysis/miscoto_bdg_458k.lp
        Running whole-community metabolic scopes...
        Community scope for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        Contributions of microbes to community scope available in output_directory/community_analysis/contributions_of_microbes.json.


        Number of metabolites producible in community: 357. 

        Reverse community scopes for all metabolic networks available in output_directory/community_analysis/rev_cscope.json and output_directory/community_analysis/rev_cscope.tsv. They higlight the producibility of metabolites by species in the community.

        --- Community scope runtime 0.72 seconds ---


        ###############################################
        #                                             #
        #    Added-value of metabolic interactions    #
        #                                             #
        ###############################################


        Added value of cooperation over individual metabolism: 23 newly reachable metabolites: 

        ...


        Added-value of cooperation written in output_directory/community_analysis/addedvalue.json

        Use the addedvalue as targets.

        Target file created with the addedvalue targets in: output_directory/community_analysis/targets.sbml
        Setting 23 compounds as targets. 


        ###############################################
        #                                             #
        #         Minimal community selection         #
        #                                             #
        ###############################################

        Running minimal community selection
        /shared/Softwares/git/miscoto/miscoto/encodings/community_soup.lp

        In the initial and minimal communities 23 targets are producible and 0 remain unproducible.

        23 producible targets:
        ...

        0 still unproducible targets:


        Minimal communities are available in output_directory/community_analysis/mincom.json 

        ######### One minimal community #########
        # One minimal community enabling the producibility of the target metabolites given as inputs
        Minimal number of bacteria in communities => 2

        ...

        ######### Key species: Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling the producibility of the target metabolites given as inputs
        Number of key species => 2

        ...

        ######### Essential symbionts: Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal communities enabling the producibility of the target metabolites given as inputs
        Number of essential symbionts => 2

        ...

        ######### Alternative symbionts: Difference between Union and Intersection #########
        # Bacteria occurring in at least one minimal community but not all minimal communities enabling the producibility of the target metabolites given as inputs
        Number of alternative symbionts => 0



        --- Mincom runtime 0.88 seconds ---

        Targets producibility are available at output_directory/producibility_targets.json
        --- Total runtime 424.18 seconds ---
        --- Logs written in output_directory/m2m_workflow.log ---


* files outputs
    * Numerous files are created in the `output_directory`, including the logs at the root of the results directory.
    
    .. code ::

        output_directory/
        ├── m2m_metadata.json
        ├── m2m_workflow.log
        ├── producibility_targets.json
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── contributions_of_microbes.json
        │   ├── mincom.json
        │   ├── rev_cscope.json
        │   ├── rev_cscope.tsv
        │   └── targets.sbml
        ├── indiv_scopes
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
        │   └── seeds_in_indiv_scopes.json
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

    ``m2m metacom`` runs the whole workflow except the reconstruction of metabolic networks. We advise to use this command to explore the metabolism of the microbial community when you already have metabolic networks.


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


More information
-----------------

Take a look at the complete tutorial in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/main/tutorials/method_tutorial>`__