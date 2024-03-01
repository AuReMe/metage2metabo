===========
m2m Outputs
===========

* m2m steps will create multiple results files in an output folder:

    .. code ::

        output_directory/
        ├── m2m_metadata.json
        ├── m2m_command.log
        ├── pgdb
        │   ├── species_1
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
        │   └── ...
        └── sbml
            ├── species_1.sbml
            └── ...
        ├── padmet
        │   ├── species_1.padmet
        │   └── ...
        └── recon_stats.tsv
        ├── indiv_scopes
        │   └── indiv_produced_seeds.json
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── contributions_of_microbes.json
        │   ├── mincom.json
        │   ├── rev_cscope.json
        │   ├── rev_cscope.tsv
        │   └── targets.sbml
        ├── producibility_targets.json

By using the ``--pwt-xml`` option, m2m will use the xml files from Pathway Tools and there will be no .dat files.

    .. code ::

        output_directory/
        ├── m2m_metadata.json
        ├── m2m_command.log
        ├── pgdb
        │   ├── species_1.xml
        │   └── ...
        └── sbml
            ├── species_1.sbml
            └── ...
        └── recon_stats.tsv
        ├── indiv_scopes
        │   └── indiv_produced_seeds.json
        │   └── indiv_scopes.json
        │   └── rev_iscope.json
        │   └── rev_iscope.tsv
        ├── community_analysis
        │   ├── addedvalue.json
        │   ├── comm_scopes.json
        │   ├── contributions_of_microbes.json
        │   ├── mincom.json
        │   ├── rev_cscope.json
        │   ├── rev_cscope.tsv
        │   └── targets.sbml
        ├── producibility_targets.json

m2m_metadata.json
-----------------

A json file summarising metadata of Metage2Metabo:

* ``m2m_args``: input arguments given to Metage2Metabo such as the command used, the path to genomes or networks, etc.

* ``tool_dependencies``: Python version, versions of Metage2Metabo and its dependencies and Clingo version. When using ``recon`` or ``workflow`` commands it will also contains Pathway Tools version.

* ``duration``: duraiton of the un of Metage2Metabo in seconds.

m2m_command.log
---------------

A log file named after the command used. For example, if you launch ``m2m metacom``, the file will be named ``m2m_metacom.log``.

This file contains all the logs return by m2m.

pgdb folder
-----------

If you use the reconstruction (with ``m2m recon`` or ``m2m workflow``), the draft metabolic network inferred by Pathway Tools for each of your species will be stored in this ``pgdb`` folder.

For each genome, if no errors occur, there will be a folder named after the genome folder name containing the Pathway-Genome Database (PGDB). The PGDB is stored using attribute-values flat files (the .dat extension file).

By using the ``--pwt-xml``, m2m will extract the xml created by MetaFlux (a module of Pathway Tools) instead of extracting the attribute-values flat .dat files.

sbml
----

After Pathway Tools metabolic network reconstruction, m2m will create SBML files from the PGDB attribute-values files. All the sbml are stored in this ``sbml`` folder.

There is one sbml for each genome given as input to m2m. These metabolic network in SBML are the input of the non-reconstruction step of m2m.

padmet
------

If you use the ``-p`` argument with the reconstruction, m2m will also create padmet file, a format used to store metabolic network information and metadata.

One padmet file is created for each genome in the PGDB folder.

recon_stats.tsv
---------------

After the reconstruction, m2m will summary the information of the draft metabolic networks in this file.

It will contain the number of genes, reactions, compounds and pathways in each metabolic networks.

indiv_scopes
------------

The indiv_scopes folder is created after the individual scopes step (in ``m2m worfklow``, ``m2m metacom`` or ``m2m iscope``). This step uses a folder containing multiples metabolic network in SBMLs and a seed file (also in SBML).

The results are stored in a json file named ``indiv_scopes.json``. The keys in this file are each metabolic network and the values are the compounds that can be produced individually by the metabolic network.

Also it can occur that seeds are produbile by individual organisms, in this case they will be listed in ``indiv_produced_seeds.json``. The keys in this file are each metabolic network and the values are the seeds that can be produced individually by the metabolic network.

The reverse of the previous iscope dictionary is stored in two files ``rev_iscope.json`` and ``rev_iscope.tsv``. The latter file is a matrix with compounds as column header and species in row. For each compounds, we have the producibility of the compounds by the species (0 not producible and 1 producible).

community_analysis
------------------

The community_analysis folder stores all the results involving the community analysis (``m2m worfklow``, ``m2m metacom``, ``m2m cscope``, ``m2m addedvalue`` or ``m2m mincom``).

comm_scopes.json
================

First step of the community analysis after the individual production analysis, the community scopes (called by ``m2m worfklow``, ``m2m metacom`` or ``m2m cscope``) shows the compounds producible by the input metabolic networks with cooperation.

The results are stored in a json with 8 keys:

* ``host_prodtargets``: if a host is given as input, contains the targets producible by the host.

* ``host_unprodtargets``: if a host is given as input, contains the targets not producible by the host.

* ``host_scope``: if a host is given as input, contains all the compounds producible by the host.

* ``com_prodtargets``: the targets producible by the community.

* ``com_unprodtargets``: the targets not producible by the community.

* ``comhost_scope``: all the compounds producible by the host + the community.

* ``com_scope``: all the compounds producible by the community.

* ``targets_producers``: for each target, the list of organisms able to produce this target. It is empty if you use ``m2m worfklow`` or ``m2m metacom`` without targets because cscope needs target to find the targets_producers.

rev_cscope.json and rev_cscope.tsv
==================================

The reverse of the ``comm_scopes.json`` scope keys are stored in two files ``rev_cscope.json`` and ``rev_cscope.tsv``. It shows for each metabolites, which species in the community can produce it.
The latter file is a matrix with compounds as column header and species in row. For each compounds, we have the producibility of the compounds by the species (0 not producible and 1 producible).

addedvalue.json
===============

After the individual scopes and the community scopes, the addedvalue (``m2m worfklow``, ``m2m metacom``, ``m2m addedvalue``), extracts the compounds that are producible by the community but not by individual organism.

The result are stored in a json file with one key ``addedvalue`` which enumerates all the compounds producible by the community but not by the individual organism.

targets.sbml
============

After the addedvalue (``m2m worfklow``, ``m2m metacom``, ``m2m addedvalue``), all the compounds that have been found by this step are stored in this sbml file. It is used as the targets file for the following step.

mincom.json
===========

Using the addedvalue or targets given by the user, the mincom step (``m2m worfklow``, ``m2m metacom`` or ``m2m mincom``) will search for the minimal community that can produce these compounds.

The results are stored in a json with 17 keys:

* ``bacteria``: organisms in the optimal solution.

* ``still_unprod``: compounds unproducible by the community.

* ``newly_prod``: compounds producible by the community.

* ``union_bacteria``: organisms from all the minimal communities.

* ``inter_bacteria``: organisms from the intersection of all the minimal communities.

* ``one_model``: results of the optimal solution.

* ``exchanged``, ``union_exchanged`` and ``inter_exchanged``: the exchanged compounds by the community, this step needs a lot of resources so it is not used in m2m. If you want to use it, use miscoto with the ``minexch`` option.

* ``key_species``: organisms from all the minimal communities.

* ``essential_symbionts``: organisms in the intersection of all the minimal communities. They are occuring in all minimal solution.

* ``alternative_symbionts``: organisms appearing in at least one minimal community but not in all.

* ``score_optimum_inter``: the optimum score found for the intersection, it corresponds to the number of organism in the minimal community.

* ``score_optimum_union``: the optimum score found for the union, it corresponds to the number of organism in the minimal community.

* ``inter_targetsproducers``: the organism that have the final reaction to produce the target in the intersection. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``union_targetsproducers``: the organism that have the final reaction to produce the target in the union. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``one_model_targetsproducers``: the organism that have the final reaction to produce the target in the optimal solution. It is a dictionary, with each target as key and the organism producing these targets as value.

contributions_of_microbes.json
==============================

A json file detailing the role of community members in the production of metabolites: which organisms have the reactions that produce the metabolites.

It contains one key per microbe in the community with several subkeys:

* ``produced_alone``: metabolties that can be produced by the organism alone.

* ``community_metabolic_gain``: metabolties that can be newly produced by the organism with the help of the community.

* ``produced_in_community``: metabolties that can be produced by the organism when it is in the community (it is composed of ``community_metabolic_gain`` and ``produced_alone``).

producibility_targets.json
--------------------------

After all these previous step, m2m (``m2m worfklow`` or ``m2m metacom``) will create this json which summarizes the producibility of each targets (either given by the user or from the addedvalue).

This json contains 12 keys:

* ``producible``: the producible compounds by the community.

* ``unproducible``: the unproducible compounds by the community.

* ``indiv_producible``: the compounds producible by individual organisms.

* ``individual_producers``: for each targets the individual organisms that can produce them.

* ``com_only_producers``: the organism that have the final reaction to produce the target but needs other organisms to produces the previous compounds needed by this final reaction. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_producible``: the compounds producible by the minimal community.

* ``key_species``: organism from all the minimal communities.

* ``essential_symbionts``: organisms in the intersection of all the minimal communities. They are occuring in all minimal solution.

* ``alternative_symbionts``: organisms appearing in at least one minimal community but not in all.

* ``mincom_inter_producers``: the organism that have the final reaction to produce the target in the intersection. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_union_producers``: the organism that have the final reaction to produce the target in the union. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_optsol_producers``: the organism that have the final reaction to produce the target in the optimal solution. It is a dictionary, with each target as key and the organism producing these targets as value.
