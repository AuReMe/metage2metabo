===========
M2M Outputs
===========

* M2M will create multiple results files in an output folder. 
    * All the possible files created by M2M are listed below:

    .. code ::

        output_directory/
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
        │   └── indiv_scopes.json
        ├── community_analysis
        │   ├── comm_scopes.json
        │   ├── addedvalue.json
        │   ├── targets.sbml
        │   ├── mincom.json
        ├── producibility_targets.json


m2m_command.log
---------------

A log file named after the command used with m2m. For example, if you launch m2m metacom, the file will be named m2m_metacom.log.

This file contains all the logs return by m2m, it summarizes the information of a run.

pgdb folder
-----------

If you use the reconstruction (with ``m2m recon`` or ``m2m workflow``), the draft metabolic network inferred by Pathway Tools for each of your species will be stored in this ``pgdb`` folder.

For each genome, if no error occurs, there will be a folder named the genome containing the Pathway-Genome Database (PGDB). The PGDB is stored using attribute-values flat files (the .dat extension file).

sbml
----

After Pathway Tools metabolic network reconstruction, m2m will create SBML from the attribute-values files (as SBML are the input of the following steps). All the sbml are stored in this ``sbml`` folder.

There is one sbml for each genome given as input to m2m.

padmet
------

If you use the ``-p`` argument with the reconstruction, m2m will also create padmet file, a format used to store metabolic network information and metadata.

recon_stats.tsv
---------------

After the reconstruction, m2m will also computes statistics about the draft metabolic networks, especially the number of genes, reactions, compounds and pathways in the metabolic network. All of these informations are stored in this file.

indiv_scopes
------------

The indiv_scopes is created after the individual scopes computation (with ``m2m worfklow``, ``m2m metacom`` or ``m2m iscope``). This computation use a folder containing multiples metabolic network in SBMLs and a seed file (also in SBML).

The results are stored in a json file ``indiv_scopes.json``. The key in this file are each metabolic network and the value are the compounds that can be produced individually by each species.

community_analysis
------------------

The community_analysis folder stores all the results involving the community analysis (``m2m worfklow``, ``m2m metacom``, ``m2m cscope``, ``m2m addedvalue`` or ``m2m mincom``).

comm_scopes.json
================

First step of the community analysis after the individual production analysis, the community scopes (called by ``m2m worfklow``, ``m2m metacom`` or ``m2m cscope``) shows the compounds producible by the input metabolic networks with cooperation.

The results are stored in a json with 8 keys:

* ``host_prodtargets``: if a host is given as input, this contains the targets producible by the host.

* ``host_unprodtargets``: if a host is given as input, this contains the targets not producible by the host.

* ``host_scope``: if a host is given as input, this contains all the compounds producible by the host.

* ``com_prodtargets``: the targets producible by the community.

* ``com_unprodtargets``: the targets not producible by the community.

* ``comhost_scope``: all the compounds producible by the host + the community.

* ``com_scope``: all the compounds producible by the community.

* ``targets_producers``: for each target, it lists the species that are involved in the last reaction to produce this target.

addedvalue.json
===============

After the individual scopes and the community scopes, the addedvalue (``m2m worfklow``, ``m2m metacom``, ``m2m addedvalue``), computes the compounds that are producible by the community but not by individual organism.

The result are stored in a json file with one key ``addedvalue`` which lists all the compounds producible by the community but not by the individual organism.

targets.sbml
============

After the addedvalue (``m2m worfklow``, ``m2m metacom``, ``m2m addedvalue``), all the compounds that have been found by this step are stored in this sbml file.

mincom.json
===========

Using the addedvalue or targets given by the user, the mincom step (``m2m worfklow``, ``m2m metacom`` or ``m2m mincom``) will search for the minimal community that can produce these compounds.

The results are stored in a json with 17 keys:

* ``bacteria``: bacteria in the optimal solution.

* ``still_unprod``: compounds unproducible by the community.

* ``newly_prod``: compounds newly producible by the community.

* ``union_bacteria``: bacteria from all the minimal communities.

* ``inter_bacteria``: bacteria from the intersection of all the minimal communities.

* ``one_model``: results of the optimal solution.

* ``exchanged``, ``union_exchanged`` and ``inter_exchanged``: the exchanged compounds by the community, this step needs a lot of resources so it is not used in m2m. If you want to use it, use miscoto with the minexch option.

* ``keystone_species``: organism from all the minimal communities.

* ``essential_symbionts``: organism in the intersection of all the minimal communities. They are occuring in all minimal solution.

* ``alternative_symbionts``: organism appearing in at least one minimal community but not in all.

* ``score_optimum_inter``: the optimum score found for the intersection, it corresponds to the number of organism in the minimal community.

* ``score_optimum_union``: the optimum score found for the union, it corresponds to the number of organism in the minimal community.

* ``inter_targetsproducers``: the organism that have the final reaction to produce the target in the intersection. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``union_targetsproducers``: the organism that have the final reaction to produce the target in the union. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``one_model_targetsproducers``: the organism that have the final reaction to produce the target in the optimal solution. It is a dictionary, with each target as key and the organism producing these targets as value.

producibility_targets.json
--------------------------

After all these previous step, m2m (``m2m worfklow`` or ``m2m metacom``) will create this json which summarizes the producibility of each targets (either given by the user or from the addedvalue).

This json contains 10 keys:

* ``producible``: the producible compounds by the community.

* ``unproducible``: the unproducible compounds by the community.

* ``indiv_producible``: the compounds producible by individual organisms.

* ``individual_producers``: for each targets the individual organisms that can produce them.

* ``com_only_producers``: the organism that have the final reaction to produce the target but needs other organisms to produces the previous compounds needed by this final reaction. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_producible``: the compounds producible by the minimal community.

* ``keystone_species``: organism from all the minimal communities.

* ``mincom_inter_producers``: the organism that have the final reaction to produce the target in the intersection. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_union_producers``: the organism that have the final reaction to produce the target in the union. It is a dictionary, with each target as key and the organism producing these targets as value.

* ``mincom_optsol_producers``: the organism that have the final reaction to produce the target in the optimal solution. It is a dictionary, with each target as key and the organism producing these targets as value.
