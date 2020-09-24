# Diabetes metagenomic study

## File contents

* `western_diet.sbml`: seeds used for m2m analyses
* `metadata.csv`: metadata file obtained from MICOM article (see `recent.csv` in [https://github.com/micom-dev/paper/tree/master/data](https://github.com/micom-dev/paper/tree/master/data))
* `addedvalue.tsv`: metabolic cooperation potential for each microbiome, resulting from m2m analyses.
* `community_scopes.tsv`: community scope for each microbiome, resulting from m2m analyses.
* `initial_community_size.tsv`: number of metabolic networks in the initial community of each sample
* `community_reduction_metrics.tsv`: number of metabolic networks in the initial community of each sample, and number of keystone species (ks), essential (es) or alternative (as) symbionts after community reduction based on the cooperation potential.
* `mds_roc_diabetes.ipynb`: jupyter notebook containing scripts used to create ROC curves and MDS figures of M2M article