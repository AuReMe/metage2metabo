# Diabetes metagenomic study

## File contents

* `western_diet.sbml`: seeds used for m2m analyses
* `metadata.csv`: metadata file obtained from MICOM article (see `recent.csv` in [https://github.com/micom-dev/paper/tree/master/data](https://github.com/micom-dev/paper/tree/master/data))
* `addedvalue.tsv`: metabolic cooperation potential for each microbiome, resulting from m2m analyses
* `community_scopes.tsv`: community scope for each microbiome, resulting from m2m analyses
* `initial_community_size.tsv`: number of metabolic networks in the initial community of each sample
* `community_reduction_metrics.tsv`: number of metabolic networks in the initial community of each sample, and number of key species (ks), essential (es) or alternative (as) symbionts after community reduction based on the cooperation potential
* `mds_roc_diabetes.ipynb`: jupyter notebook containing scripts used to create ROC curves and MDS figures of M2M article
* `initial_community_composition.tsv`: phylum-level composition of initial communities in individuals
* `community_reduction_composition.tsv`: phylum-level composition of key species, alternative and essential symbionts after community reduction optimising the producibility of the cooperation potential (added_value) or the entire set of producible compounds in the community (full_targets)
* `t2d_m2m_size.R`: R script to generate figures comparing various metrics between statuses, either for MHD, SWE or all samples (change the focus variable)
* `t2d_m2m_target_producers.tsv`: analyse butyrate target producers and produce pie plots showing taxonomic diversity of the full metagenomic dataset
* `t2d_m2m_tax.R`: analyse taxonomic composition of different groups (key, essential etc). used in fig 4 and supp figs
* `abundance.mat`: abundance matrix
* `individual_scopes.tsv`: individual producers of all metabolites
* `t2d_m2m`: R script to perform analyses on community scopes and cooperation potential compositions and generate figures of the paper
* `community_scopes_long.tsv`: long version of community scopes for each sample
* `community_scopes_by_family.tsv`: community scopes for each sample, summarised by families of metabolites
* `addedvalue_long.tsv`: long version of cooperation potential for each sample
* `ko_by_genome.tsv`: matrix describing the KO content of each genome according to the Eggnog-mapper annotation
* `t2d_m2m_div.R`: R script that analyses the diversity/richness of samples as well as their metabolic distance