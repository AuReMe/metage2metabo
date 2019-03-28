===
m2m
===

Metage2metabo is a Python3 tool to perform graph-based metabolic analysis starting from annotated genomes (**reference genomes or metagenome-assembled genomes**.
It uses **Pathway Tools** in a automatic and parallel way to *reconstruct metabolic networks* for a large number of genomes.
The obtained metabolic networks are then *analyzed individually and collectively* in order to get the *added value of metabolic cooperation in microbiota* over individual metabolism, and to *identify and screen interesting organisms* among all.

m2m can be used as a whole workflow (``m2m workflow``) or steps can be performed individually (``m2m recon``, ``m2m iscope``, ``m2m cscope``, ``m2m addedvalue``, ``m2m mincom``, ``m2m seeds``).

This project is licensed under the GNU General Public License - see the `LICENSE.md <https://github.com/AuReMe/metage2metabo/blob/master/LICENSE>`__ file for details.

m2m's documentation
===================

.. toctree::
    :maxdepth: 2
    :glob:

    install
    command
    tutorial
    suppinformation
    api


Additional features
===================

M2M relies on packages that can also be used independantly with more features:

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
