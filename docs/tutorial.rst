============
m2m Tutorial
============
Test data is avaible in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/test>`__.
It contains enough data to run the different subcommands.

m2m recon
---------
``m2m recon`` runs metabolic network reconstruction for all annotated genomes, using Pathway Tools. It can be done with multiple CPUs, in which case the number of allocated/available CPU has to be given as a optional argument.

It uses the following mandatory inputs (run ``m2m recon --help`` for optional arguments):

-g directory           directory of annotated genomes
-s file                seeds SBML file
-o directory           output directory for results

Optional arguments:

-c int           number of CPU for multi-processing
--clean          option to rerun every reconstruction 
                 even if found in ptools-local

The inputs genomic data has to follow a strict strure, that can be observed in the `recon_data` directory of tests in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/test>`__. It is reproduced below:

::

    recon_data
        ├── organism_1          
        │   └── organism_1.gbk
        ├── organism_2          
        │   └── organism_2.gbk
        ├── organism_3          
        │   └── organism_3.gbk
        ..
        └── organism_n         
            └── organism_n.gbk

.. code:: sh

    m2m recon -g test/recon_data -o output_directory -c cpu_number [--clean]

* standard output
    .. code:: 

        ######### Running metabolic network reconstruction with Pathway Tools #########
        tca_cycle_ecolicyc (at xxxx/ptools-local/pgdbs/user/tca_cycle_ecolicyc) has been removed.
        425hofleacyc (at xxxx/ptools-local/pgdbs/user/425hofleacyc) has been removed.
        fatty_acid_beta_oxydation_icyc (at xxxx/ptools-local/pgdbs/user/fatty_acid_beta_oxydation_icyc) has been removed.
        404rhizobiumcyc (at xxxx/ptools-local/pgdbs/user/404rhizobiumcyc) has been removed.
        ~~~~~~~~~~Creation of input data from Genbank/GFF~~~~~~~~~~
        Checking inputs for tca_cycle_ecoli: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for tca_cycle_ecoli.
        Checking inputs for fatty_acid_beta_oxydation_I: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for fatty_acid_beta_oxydation_I.
        ~~~~~~~~~~Inference on the data~~~~~~~~~~
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho recon_data//tca_cycle_ecoli/
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho recon_data//fatty_acid_beta_oxydation_I/
        ~~~~~~~~~~Check inference~~~~~~~~~~

        2 builds have passed!

        ~~~~~~~~~~Creation of the .dat files~~~~~~~~~~
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load recon_data//tca_cycle_ecoli//dat_creation.lisp
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load recon_data//fatty_acid_beta_oxydation_I//dat_creation.lisp
        ~~~~~~~~~~Check .dat ~~~~~~~~~~
        tca_cycle_ecolicyc: 23 on 23 dat files create.
        fatty_acid_beta_oxydation_icyc: 23 on 23 dat files create.
        ~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~
        ~~~~~~~~~~Moving result files~~~~~~~~~~
        ~~~~~~~~~~The script have finished! Thank you for using it.
        ######### Creating SBML files #########
        PGDB created in output_directory//pgdb
        SBML files created in output_directory//sbml

        Here the ``--clean`` option was used. The output shows that PGDB are created with Pathway Tools. Then the .dat files are extracted and used to build SBML files of the metabolic models. 
* files outputs
    * In `output_directory/pgdb` are found the .dat files of Pathway Tools. The corresponding SBMLs are in `output_directory/sbml`. The structure of the output directory after this ``recon`` command is shown below :

    ::

        output_directory/
        ├── pgdb
        │   ├── fatty_acid_beta_oxydation_I
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
        │   └── tca_cycle_ecoli
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
        └── sbml
            ├── fatty_acid_beta_oxydation_I.sbml
            └── tca_cycle_ecoli.sbml

        * Finally, in the input directory, some files are also generated automatically by Pathway Tools
        ::
            
            recon_data/
            ├── fatty_acid_beta_oxydation_I
            │   ├── dat_creation.lisp
            │   ├── fatty_acid_beta_oxydation_I.gbk
            │   ├── genetic-elements.dat
            │   ├── organism-params.dat
            │   └── pathologic.log
            └── tca_cycle_ecoli
                ├── dat_creation.lisp
                ├── genetic-elements.dat
                ├── organism-params.dat
                ├── pathologic.log
                └── tca_cycle_ecoli.gbk


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

It uses the following mandatory inputs (run ``m2m mincom --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file                targets SBML file
-o directory           output directory for results

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
    * In `output_directory/indiv_scopes/indiv_scopes.json`. A json file that can be easily loaded as a dictionary (or humanly read as it it) that contains the set of reachable metabolites for each organism. /!\\ Warning: the seeds are included in the scopes, hence they will never be empty. 

cscope
*******

It uses the following mandatory inputs (run ``m2m mincom --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file                targets SBML file
-o directory           output directory for results

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
    * In addition, a json file with the results is created in `output_directory/community_analysis/indiv_scopes.json`.

addedvalue
**********

``m2m addedvalue`` uses the previously two subcommands to compute the added value of combining metabolisms in the microbiota (i.e. consider metabolic cooperation) with respect to studying individually the metabolism of each organism. 
It uses the following mandatory inputs (run ``m2m addedvalue --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-o directory           output directory for results

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
    * A targets SBML file is generated. It can be used with `` m2m mincom`` . The json files associated to ``iscope`` and ``cscope`` are also produced.

    ::

        output_directory/
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── miscoto_om6hubmz.lp
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

.. code:: sh

    m2m mincom -n toy_bact -s metabolic_data/seeds_toy.sbml -t metabolic_data/targets_toy.sbml -o output_directory/

* standard output
    .. code::

        ######### Creating metabolic instance for the whole community #########
        Created instance in output_directory/community_analysis/miscoto_36t8lqe_.lp
        Running minimal community selection
        Community scopes for all metabolic networks available in output_directory/community_analysis/comm_scopes.json
        ######### One minimal community #########
        # One minimal community enabling to produce the target metabolites given as inputs
        Minimal number of bacteria in communities = 13
        GCA_003437905
        GCA_003437255
        GCA_003437055
        GCA_003437885
        GCA_003437815
        GCA_003437595
        GCA_003437375
        GCA_003438055
        GCA_003437665
        GCA_003437945
        GCA_003437295
        GCA_003437195
        GCA_003437715
        ######### Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling to produce the target metabolites given as inputs
        Union of bacteria in minimal communities = 17
        GCA_003437345
        GCA_003437905
        GCA_003437255
        GCA_003437055
        GCA_003437175
        GCA_003437885
        GCA_003437325
        GCA_003437815
        GCA_003437595
        GCA_003437375
        GCA_003438055
        GCA_003437665
        GCA_003437945
        GCA_003437295
        GCA_003437195
        GCA_003437715
        GCA_003437785
        ######### Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal community enabling to produce the target metabolites given as inputs
        Intersection of bacteria in minimal communities = 12
        GCA_003437905
        GCA_003437255
        GCA_003437055
        GCA_003437885
        GCA_003437815
        GCA_003437595
        GCA_003437375
        GCA_003438055
        GCA_003437665
        GCA_003437295
        GCA_003437195
        GCA_003437715

    This output gives the result of minimal community selection. It means that for producing the 119 metabolic targets, a minimum of 13 bacteria out of the 17 is required. One example of such minimal community is given. In addition, the whole space of solution is studied. All bacteria (17) occur in at least one minimal community. Finally, the intersection gives the following information: a set of 12 bacteria occurs in each minimal communtity. This means that these 12 bacteria are needed in any case, and that any of the remaining 5 bacteria can complete the missing function(s).
* files outputs
    * As for other commands, a json file with the results is produced in ``output_directory/community_analysis/comm_scopes.json``

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

You can run the workflow analysis with the two genbanks files available in the `Github repository <https://github.com/AuReMe/metage2metabo/tree/master/test>`__ (`workflow_data`). Two genomes are available in the compressed archive workflow_genomes.tar.gz. The archive has to be uncompressed before testing.

    .. code:: sh

        m2m workflow -g workflow_genomes -s workflow_data/seeds_workflow.sbml -o output_directory/

* standard outputs
    .. code ::

        ######### Running metabolic network reconstruction with Pathway Tools #########
        ~~~~~~~~~~Creation of input data from Genbank/GFF~~~~~~~~~~
        Checking inputs for GCA_900318805: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for GCA_900318805.
        Checking inputs for GCA_900315385: missing organism-params.dat; genetic-elements.dat; dat_creation.lisp. Inputs file created for GCA_900315385.
        ~~~~~~~~~~Inference on the data~~~~~~~~~~
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_900318805/
        pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho workflow_genomes/GCA_900315385/
        ~~~~~~~~~~Check inference~~~~~~~~~~

        2 builds have passed!

        ~~~~~~~~~~Creation of the .dat files~~~~~~~~~~
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load workflow_genomes/GCA_900318805//dat_creation.lisp
        pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load workflow_genomes/GCA_900315385//dat_creation.lisp
        ~~~~~~~~~~Check .dat ~~~~~~~~~~
        gca_900318805cyc: 23 on 23 dat files create.
        gca_900315385cyc: 23 on 23 dat files create.
        ~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~
        ~~~~~~~~~~Moving result files~~~~~~~~~~
        ~~~~~~~~~~The script have finished! Thank you for using it.
        ######### Creating SBML files #########
        ######### Running individual metabolic scopes #########
        Individual scopes for all metabolic networks available in output_directory//indiv_scopes/indiv_scopes.json
        2 metabolic models considered.
        29 metabolites in core reachable by all organisms (intersection)
        37 metabolites reachable by individual organisms altogether (union), among which 26 seeds (growth medium)
        max metabolites in scope 36
        min metabolites in scope 30
        average number of metabolites in scope 33.00 (±4.24)
        ######### Creating metabolic instance for the whole community #########
        Created instance in output_directory/community_analysis/miscoto_ena_9l33.lp
        Running whole-community metabolic scopes
        Community scopes for all metabolic networks available in output_directory//community_analysis/comm_scopes.json
        Added value of cooperation over individual metabolism: 25 newly reachable metabolites:
        M_2__45__PG_c, M_METHYL__45__GLYOXAL_c, M_D__45__SEDOHEPTULOSE__45__7__45__P_c, M_NITRITE_c, M_DIHYDROXY__45__ACETONE__45__PHOSPHATE_c, M_FRUCTOSE__45__16__45__DIPHOSPHATE_c, M_GAP_c, M_RIBOSE__45__5P_c, M_RIBULOSE__45__5P_c, M_CPD__45__12079_c, M_G3P_c, M_PHOSPHORIBOSYL__45__FORMIMINO__45__AICAR__45__P_c, M_NADH_c, M_PRPP_c, M_DPG_c, M_3__45__P__45__HYDROXYPYRUVATE_c, M_PHOSPHORIBULOSYL__45__FORMIMINO__45__AICAR__45__P_c, M_PHOSPHORIBOSYL__45__AMP_c, M_L__45__LACTATE_c, M_ERYTHROSE__45__4P_c, M_PHOSPHORIBOSYL__45__ATP_c, M_D__45__LACTATE_c, M_XYLULOSE__45__5__45__PHOSPHATE_c, M_BETA__45__D__45__FRUCTOSE_c, M_FRUCTOSE__45__6P_c
        Setting these 25 as targets
        Running minimal community selection
        Community scopes for all metabolic networks available in output_directory//community_analysis/comm_scopes.json
        ######### One minimal community #########
        # One minimal community enabling to produce the target metabolites given as inputs
        Minimal number of bacteria in communities = 2
        GCA_900318805
        GCA_900315385
        ######### Union of minimal communities #########
        # Bacteria occurring in at least one minimal community enabling to produce the target metabolites given as inputs
        Union of bacteria in minimal communities = 2
        GCA_900318805
        GCA_900315385
        ######### Intersection of minimal communities #########
        # Bacteria occurring in ALL minimal community enabling to produce the target metabolites given as inputs
        Intersection of bacteria in minimal communities = 2
        GCA_900318805
        GCA_900315385
* files outputs
    * Numerous files are created in the output_directory
    
    .. code ::

        output_directory/
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── mincom.json
        │   ├── miscoto_ena_9l33.lp
        │   ├── miscoto_ena_9l33__tgts.lp
        ├── indiv_scopes
        │   └── indiv_scopes.json
        ├── pgdb
        │   ├── GCA_900315385
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
        │   └── GCA_900318805
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
        └── sbml
            ├── GCA_900315385.sbml
            └── GCA_900318805.sbml


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
