# Changelog

# Metage2Metabo v1.6.0 (2024-03-01)

WARNING: change for individual and community scopes:
* for individual scope: only seeds that are producible, ie associated with activated reactions are now shown in the results. Additional information on producible, non-producible and absent seeds are available in the json output. This requires MeneTools version >= `3.4.0`.
* for community scope: use of `miscoto focus` for computation in order to retrieve the metabolite producers in community. This requires MiSCoTo version >= `3.2.0`.

## Add

* Test in m2m_analysis to check if combination of powergraph predicted by bubbletools is the same as the ones found in minimal communities.
* Creation of boolean equation summarizing the powergraph if it is simple enough.
* New option `--target-com-scope` to use all the community scope as targets for minimal community prediction (issue #21).
* Add a function to modify xml created by Pathway Tools (issue #60).
* Troubleshooting page in readthedocs (issue #24 and #25).
* Metadata json file created when using command line (both for `m2m` and `m2m_analysis`).
* Essential and alternative symbiont in `producibility_targets.json` (issue #19).
* Test for host for m2m metacom.

## Fix

* Fix potential error when reading taxon_id file (issue #22).
* Fix issue with etree in reconstruction.
* Issues with readthedocs.
* GitHub Actions not working.

## Modify

* Better deal with seeds that are absent from networks or not produced in iscope and also the seeds that are produced through interactions in cscope (issue #53).
* Exit m2m_analysis when there are unproducible targets (issue #23).
* Check forbidden characters during targets file creation.
* Do not allow abbreviation in command arguments by argparse (issue #54).
* Sanitize use of tarfile extractall (issue #55).
* Remove unused dependency.
* First step in replacing `pkg_resources` (which will become deprecated in the future) with `importlib.metadata` or import of \_\_version\_\_.
* Move from `setup.py`/`setup.cfg` to `pyproject.toml`.
* Update tutorial.
* Update docs and readme.
* Update license year and affiliation.
* Move support for CI to Python 3.8 and 3.9.
* Remove import in `__init__.py` to avoid loading m2m_analysis dependencies when using only m2m.
* Avoid using os.unlink() function on Windows as it can lead to Permission Error.

# Metage2Metabo v1.5.4 (2023-05-25)

## Fix

* Issues in html powergraph creation:
  * some nodes are not indicated as alternative symbionts.
  * colors are not applied to some nodes.
  * the names of the nodes contain unicode ints.

## Modify

* Update license year.

# Metage2Metabo v1.5.3 (2022-09-26)

## Fix

* Correctly fix the issue with the colors used for the powergraph.

# Metage2Metabo v1.5.2 (2022-09-23)

## Modify

* Change the shapes of the nodes in the html output of `m2m_analysis`: circle for essential symbionts and square for alternative symbionts.

## Fix

* Issue with the number of colors used to color the taxon in the powergraph.

# Metage2Metabo v1.5.1 (2022-09-21)

## Fix

* issue in m2m_analysis where clyngor module uses the clingo module producing [ValueError](https://github.com/Aluriak/PowerGrASP/issues/1).
* issue in m2m_analysis when incorrect paths for seed and/or target files were given as input. 

# Metage2Metabo v1.5.0 (2021-03-17)

This release focuses on `m2m_analysis` by adding new option, new output files,more documentation and refactoring some functions and output files.

There is also a complete refactoring on the structure of m2m repository: the two big scripts (`m2m_workflow.py` and `m2m_analysis.py`) are now split in multiple scripts in two folders ( `m2m` and `m2m_analysis`).

## Add

* New option for `m2m_analysis`: `--level` to select the taxonomy level (phylum, class, order, family, genus or species) that will be used in the following steps of m2m_analysis, by default, it is phylum.
* HTML output for `m2m_analysis`. HTML visualization is interactive and nodes are coloured (by key species types or taxon).
* Error message with incorrect SBML file format.
* Log files for `m2m_analysis`.
* New tests for `m2m_analysis`.
* Documentation about `m2m_analysis` output files.
* Citation to eLife article.
* Tutorial with a jupyter notebook about `m2m metacom` methods.

## Modify

* Add colours (by key species types or taxon) to `m2m_analysis` svg.
* Modify `m2m_analysis` output files.
* Replace GPL license by LGPL license.
* Complete refactoring of m2m repository structure: `m2m_workflow.py` and `m2m_analysis.py` have been split into multiple scripts contained inside `m2m` and `m2m_analysis` folders.

## Remove

* Delete `m2m_analysis stats` because it is no more useful with the changes made to `m2m_analysis`.

# Metage2Metabo v1.4.1 (2020-12-17)

## Add

* Clarified logs when some target metabolites are also part of the seeds. Indications to which file one must refer in order to check the individual producibility of such seeds by organisms.

# Metage2Metabo v1.4.0 (2020-12-14)

**Warning**: the version 1.4.0 of m2m needs updated version of its dependencies:
- menetools >= 3.1.0.
- miscoto >= 3.0.3.
- mpwt >= 0.6.0. 

## Add

* Multiprocessing for `m2m iscope` and `m2m metacom`, wich can be called using the options `-c`. 
* New option `--pwt-xml` for `m2m recon` and `m2m workflow`: this option extracts the XML file created by MetaFlux (a module of Pathway Tools) so it can be used as input for m2m instead of extracting the `.dat` files and using padmet to create the SBML input (using mpwt 0.6.0).
* Scripts and data for diversity analysis in article data.
* New output file "indiv_produced_seeds.json" showing if seeds are producible by individual organisms (using menetools 3.1.0).

## Fix

* Issue with seeds in individual_producers.
* SBML and padmet stats computation for `m2m recon` and `m2m workflow` so it can count correctly reactions, metabolites and genes in XML from MetaFlux.

## Modify

* Update article data, scripts and jupyter notebook.
* Update readme and docs.
* Update requirements.txt and delete redundant files.
* Replace `keystone species` by `key species` (also modified in miscoto 3.0.3).

# Metage2Metabo v1.3.5 (2020-11-18)

## Add

* Windows and MacOS tests in GitHub Actions.
* Windows compatibility for `m2m metacom` using menetools 3.0.2 and miscoto 3.0.2.

# Metage2Metabo v1.3.4 (2020-10-30)

## Fix

* Issue with warning message (about compounds being both  in seeds and targets) always showing.

# Metage2Metabo v1.3.3 (2020-10-06)

## Fix

* Seeds file validation broke for m2m options that did not require seeds (no "seeds" in args namespace", e.g `m2m recon`)

## Add

* data and scripts associated to the manuscript in the Github repository

# Metage2Metabo v1.3.2 (2020-09-18)

## New features

* Provide a proxy to estimate functional redundancy in metabolite producibility by displaying, for each metabolite, all the metabolic networks that can produce it in the initial community. Outputs are a matrix in `indiv_scopes/rev_iscope.tsv` and in `indiv_scopes/rev_iscope.json`
* Provide a warning message + explanations in logs when one or several targets are also seeds. 

## Fix

* Enhance the readability of logs to clarify the producibility results of targets when a target file is given in inputs.
* Show all producible metabolites in `community_analysis/mincom.json` dict key "producible" so that when a target is also a seed, it appears under that key together with all other producible targets (a seed that is a target does not appear under the key `newly_producible`). Associated to updates in MiSCoTo (v2.1.2: this version is now a requirement).
* Fix the logs for the use case above
* Singularity image in Singularity Hub for `m2m` commands that do not require Pathway Tools
* Remove host keys in `community_analysis/comm_scopes.json` if no host is provided. 

## Doc

* New outputs in `indiv_scopes` to estimate functional redundancy
* Explanations about the new singularity image in Singularity Hub

## CI 

* Test for metabolite producers (reverse iscopes)
* Tests to account for the presence of targets that are also seeds

# Metage2Metabo v1.3.1 (2020-07-28)

## Add

* Targets option (`-t`)  for cscope.
* New output file called `producibility_targets.json`, showing the producibility of each targets and by which organisms.
* Final producers for the target using miscoto 2.1.1. The final producers are the organisms able to produce the target in full community or in minimal communities.
* New page in the readthedocs with all the output of m2m.
* Seeds from the diabetes metagenomic analysis of the paper.

## Modification

* Remove sbml conversion, as it seems that menetools and miscoto can handle SBML 3.

## Fix

* Issue with m2m test and `-t` option.

## Update

* Doc and README and fix numerous typos.

## CI

* Test for cscope, the new output file and final producers given by miscoto.

# Metage2Metabo v1.3.0 (2020-07-07)

## Add

* Targets option (`-t`) for metacom and workflow. By using this option, m2m will replace the added-value targets by the user targets.
* Use Miscoto function to_json to create mincom json, this fix issue with the formatting of the json. M2M now requires Miscoto >= 2.1.0.

## Fix

* Oog.jar is now also check in `__main_analysis__` to show if there is an error with the Oog.jar file earlier.
* Issue with logs especially mpwt logs not showing off.

## Update

* Doc and README and fix numerous typos.

## CI
* Do not launch the CI if the commit modify only the readme or the docs.

# Metage2Metabo v1.2.1 (2020-06-18)

## Add

* Error messages during compounds sbml creation.

## Fix

* Issue with mpwt results check.

# Metage2Metabo v1.2.0 (2020-04-22)

## New features

* `m2m metacom` to run the whole workflow *except* metabolic network reconstruction
* Add a log file by default in the results directory
* Add a quiet option `- q` that write in the log file and streams to stdout logs from warning level
* Create a json output file when calling `m2m addedvalue`
* Add a new option `-p` in recon to create padmet files. 
* Create a `targets.sbml` file when running `addedvalue` in  `output folder/community_analysis`

## Fix

* Typos in doc, README
* Logger issues
* Build of the package in ReadTheDoc
* Behaviour of the programme when using empty metabolic networks
* Retrieve the caught exceptions from the dependencies to give more information for troubleshooting (especially useful for mpwt and Pathway Tools)

## Update

* Doc and README with 
    * more details on the modelling of producibility
    * documentation of new functionalities (`metacom`, quiet option...)
* Compatibility with mpwt 0.5.4 


## CI

* Use of Github Actions instead of Travis
* Test for `m2m metacom`
* Fix and customise existing tests

# Metage2Metabo v1.1.6 (2020-01-13)

## Fix

* compatibility with mpwt 0.5.3.
* typos.
  
# Metage2Metabo v1.1.5 (2019-11-19)

## Add

* check for powergrasp and ete3 installations in m2m_analysis.

## Fix

* compatibility with padmet 4.0.0.

# Metage2Metabo v1.1.4 (2019-10-29)

## Fix

* issues with malformatted SBML (issue #7).

# Metage2Metabo v1.1.3 (2019-10-28)

## Fix

* m2m tries to create a target sbml from an empty added value.

# Metage2Metabo v1.1.2 (2019-10-23)

## Add

* test for sbml conversion.
* tests for m2m_analysis.
* article seeds in article_data.

## Fix

* issue with m2m recon and workflow when reconstructing a metabolic network from only one organism.

## Update

* requirements.
* Singularity and Docker recipes.
* travis yml file with m2m_analysis dependencies.

# Metage2Metabo v1.1.1 (2019-10-23)

## Add

* log to m2m_analysis.
* keystone species to m2m log.

## Fix

* m2m_analysis readthedocs.
* numerous issues with m2m_analysis.

# Metage2Metabo v1.1.0 (2019-10-23)

## Add

* **m2m_analysis**: optimal solutions enumeration, graph representation and compression into a powergraph.
* padmet creation after Pathway Tools reconstruction (wtih -p argument for m2m recon and workflow).
* statistics on genes/reactions/compounds in GSM after the reconstruction.

## Fix

* miscoto and clyngor memory issue with solutions enumeration.

## Update

* readthedocs by adding m2m_analysis.

# Metage2Metabo v1.0.0 (2019-10-23)

## Update/Fix

* replace pyasp with [clyngor](https://github.com/Aluriak/clyngor).

# Metage2Metabo v0.0.10 (2019-06-03)

First release, uses Pyasp as a support for ASP computing.
