Metage2Metabo (M2M)
===================

Metage2metabo is a Python3 tool to perform graph-based metabolic analysis starting from annotated genomes (**reference genomes or metagenome-assembled genomes**.
It uses **Pathway Tools** in a automatic and parallel way to *reconstruct metabolic networks* for a large number of genomes.
The obtained metabolic networks are then *analyzed individually and collectively* in order to get the *added value of metabolic cooperation in microbiota* over individual metabolism, and to *identify and screen interesting organisms* among all.

M2M can be used as a whole workflow (``m2m workflow``) or steps can be performed individually (``m2m recon``, ``m2m iscope``, ``m2m cscope``, ``m2m addedvalue``, ``m2m mincom``, ``m2m seeds``).

This project is licensed under the GNU General Public License - see the `LICENSE.md <https://github.com/AuReMe/metage2metabo/blob/master/LICENSE>`__ file for details.

General information about the modelling
========================================

M2M has two main dependencies for modelling metabolic networks: `MeneTools <https://github.com/cfrioux/MeneTools>`__ and `Miscoto <https://github.com/cfrioux/miscoto>`__. Accordingly metabolic models in M2M follow the producibility in metabolic networks as defined by the `network expansion <http://www.ncbi.nlm.nih.gov/pubmed/15712108>`__ algorithm.
Mainly, two rules are followed:
* a *recursive rule*: the products of a reactions are producible if **all** reactants of this reaction are themselves producible
* an *initiation rule*: producibility is initiated by the presence of nutrients, called *seeds*. 

A metabolite that is producible from a set of nutrients is described as being "in the scope of the seeds".
The computation is made using logic solvers (Answer Set Programming). The present modelling ignores the stoichiometry of reactions (2A + B --> C is considered equivalent to A + B --> C), and is therefore suited to non-curated or draft metabolic networks, as the ones built using M2M with the PathoLogic software of `Pathway Tools <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5036846/pdf/bbv079.pdf>`__ handled by `Mpwt <https://github.com/AuReMe/mpwt>`__. Many works have relied on network expansion to study organisms (`here <http://doi.wiley.com/10.1111/tpj.12627>`__, `here <https://dx.plos.org/10.1371/journal.pcbi.1000049>`__ or `there <http://dx.plos.org/10.1371/journal.pcbi.1005276>`__) and communities (`here <https://academic.oup.com/bioinformatics/article/34/17/i934/5093211>`__, `here <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4786-7>`__, or `here <https://www.ncbi.nlm.nih.gov/pubmed/18546499>`__). It has been `compared <http://www.ncbi.nlm.nih.gov/pubmed/19425125>`__, `combined <https://www.cambridge.org/core/product/identifier/S1471068418000455/type/journal_article>`__ to steady-state modelling (Flux Balance Analysis). 

Documentation of M2M
====================

.. toctree::
    :maxdepth: 2
    :glob:

    install
    command
    tutorial
    suppinformation
    api
    m2m_analysis


Additional features
===================

M2M relies on packages that can also be used independantly with more features and additional options:

* `mpwt <https://github.com/AuReMe/mpwt>`__: command-line and multi-process solutions to run Pathway Tools. Suitable to multiple reconstruction, for example genomes of a microbiota

* `menetools <https://github.com/cfrioux/MeneTools>`__: individual metabolic capabilities analysis using graph-based producibility criteria

* `miscoto <https://github.com/cfrioux/miscoto>`__: community selection and metabolic screening in large-scal microbiotas, with or without taking a host into account

Authors
=======

`Clémence Frioux <https://cfrioux.github.io/>`__ and `Arnaud Belcour <https://arnaudbelcour.github.io/blog/>`__, Univ Rennes, Inria, CNRS, IRISA, Rennes, France.

Acknowledgement
===============

People of Pathway Tools (SRI International) for their help integrating Pathway Tools with command line and multi process in the `mpwt <https://github.com/AuReMe/mpwt>`__ package, used in M2M.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
