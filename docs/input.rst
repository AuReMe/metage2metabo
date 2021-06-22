=========
m2m input
=========

Genomes or metabolic networks
-----------------------------

If you have annotated genomes and no metabolic networks, you can reconstruct draft metabolic networks using ``m2m`` with the ``m2m workflow`` or the ``m2m recon`` commands.
The genomes are expected to be in GenBank format and the input is the directory containing subdirectories with the corresponding files.
In addition, the name of the directory must match the name of the genbank file.
Example:

::

input_folder
├── organism_1
│ └── organism_1.gbk
├── organism_2
│ └── organism_2.gbk
├── organism_3
│ └── organism_3.gbk
..
└── organism_n
└── organism_n.gbk


If you already have metabolic networks, you can analyse them with the command ``m2m metacom`` that will perform the whole workflow except for the metabolic reconstruction part.
The metabolic networks must be in SBML format. The input is a folder containing multiples SBML files.
Example:

::

input_folder
├── organism_1.sbml
├── organism_2.sbml
├── organism_3.sbml
..
└── organism_n.sbml

Seeds
-----

Seeds are a set of metabolites that define the nutritional conditions of the community.
These metabolites can be components of the growth medium or co-factors of complex cycles (known as currency metabolites).

For example, this is a list of common currency metabolites (from Kim et al. 2015):

- proton
- water
- oxygen molecule
- NADP+
- NADPH
- ATP
- diphosphate
- carbon dioxide
- phosphate
- ADP
- coA
- UDP
- NAD+
- NADH
- AMP
- ammonia
- hydrogen peroxide
- oxidized electron acceptor
- reduced electron acceptor
- 3-5-ADP
- GDP
- carbon monoxide
- GTP
- FAD

The seed metabolites can be found in the literature.
They can be also found from metabolic databases, for example metaboltites from `nutrition in the VMH <https://www.vmh.life/#nutrition>`__ or `growth media from MetaCyc <https://metacyc.org/META/new-image?object=Growth-Media>`__ (in this latter case, one should used the metabolites in the ``Constituents`` column).

Once a list of metabolites has been designed, these metabolites must be converted to a list of IDs of the metabolic database corresponding to your metabolic networks.
For example, in the VMH seeds ethanol is ``etoh``. In the MetaCyc database, the ID of ``ethanol`` is ``ETOH``.
Then the ID must be checked with the ID in the SBML files of the metabolic networks. In this example ``ETOH`` is associated to ``M_ETOH_c`` in the SBML file in the species field (``<species id="M_ETOH_c" name="ETOH" compartment="c"/>``).
The ``M_`` corresponds to metabolite (another possibility for this prefix is the ``R_`` for reaction). And ``_c`` corresponds to the cytosol compartment.

Targets (Optional)
------------------

Targets are a set of metabolites whose producibility by the community will be checked and for which minimal communities will be computed.
It is possible to give targets to Metage2Metabo with any of the commands ``m2m cscope``, ``m2m mincom`` and ``m2m metacom``.
Adding targets will replace the metabolites from the cooperation potential (the ``addedvalue``) as a metabolic objective whether these metabolites are producible either by individual organisms (with ``iscope``), by the community (with ``cscope``) and compute minimal communities producing these metabolites (with ``mincom``).

The same conversion method than the one for the seeds is needed to create a useful targets file.
In addition, one can customise the set of predefined targets from the cooperation potential by modifying the file ``community_analysis/targets.sbml`` and setting this file as targets with the ``-t`` argument.

Bibliography
------------

Kim, T., Dreher, K., Nilo-Poyanco, R., Lee, I., Fiehn, O., Lange, B. M., Nikolau, B. J., Sumner, L., Welti, R., Wurtele, E. S., & Rhee, S. Y. (2015). Patterns of metabolite changes identified from large-scale gene perturbations in Arabidopsis using a genome-scale metabolic network. Plant physiology, 167(4), 1685–1698. https://doi.org/10.1104/pp.114.252361