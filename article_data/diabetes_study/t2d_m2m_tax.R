# load packages
library(reshape2)
library(dplyr)
library(ggplot2)
library(ape)
library(vegan)
library(viridis)
library(RColorBrewer)

######### CONFIG ###########
focus = "mhd" # can be "mhd" or "swe" or "all"
############################

metadata_file = 'metadata.csv'
tax_file = 'mgs_specI.tax'
taxcom_sample_ks_file = 'community_reduction_composition.tsv'
taxcom_sample_init_file = 'initial_community_composition.tsv'
metrics_file_av = 'community_reduction_metrics.tsv'

res_dir = paste0('figures_results_', focus, "_tax", "/")
dir.create(file.path(res_dir), showWarnings = FALSE)

# we want to drop samples that had less than 1000 reads passing quality checks
samples_to_be_dropped = c("NG5425286", "NG5425485", "MH0035", "MH0037", "MH0038", "MH0047", "MH0054", "MH0057", "MH0058", "MH0064", "MH0065", "MH0067", "MH0073", "MH0077", "MH0082", "MH0084")

# read input files
metrics_adval = as.data.frame(read.table(metrics_file_av, header=TRUE, sep="\t"))
taxcom_sample_init = as.data.frame(read.table(taxcom_sample_init_file, header=TRUE, sep="\t"))
taxcom_sample_ks = as.data.frame(read.table(taxcom_sample_ks_file, header=TRUE, sep="\t"))

# read metadata
metadata = as.data.frame(read.table(metadata_file, header=TRUE, sep=","))
tax = as.data.frame(read.table(tax_file, header=FALSE, sep="\t"))
names(tax) = c("genome","K","P","C","O","F","G","S")

# remove samples that we dropped in each file
metadata = metadata[! metadata$sample %in% samples_to_be_dropped, ]
taxcom_sample_init = taxcom_sample_init[! taxcom_sample_init$sample %in% samples_to_be_dropped, ]
taxcom_sample_ks = taxcom_sample_ks[! taxcom_sample_ks$sample %in% samples_to_be_dropped, ]
metrics_adval = metrics_adval[! metrics_adval$sample %in% samples_to_be_dropped,]


mhd_samples = metadata[metadata$subset == 'MHD', 'sample']
swe_samples = metadata[metadata$subset == 'SWE', 'sample']

if (focus=="mhd"){
    # drop the SWE samples
    metadata = metadata[! metadata$sample %in% swe_samples,]
    taxcom_sample_init = taxcom_sample_init[! taxcom_sample_init$sample %in% swe_samples,]
    taxcom_sample_ks = taxcom_sample_ks[! taxcom_sample_ks$sample %in% swe_samples,]
    metrics_adval = metrics_adval[! metrics_adval$sample %in% swe_samples,]
} else if (focus=="swe"){
    # drop the MHD samples
    metadata = metadata[! metadata$sample %in% mhd_samples,]
    taxcom_sample_init = taxcom_sample_init[! taxcom_sample_init$sample %in% mhd_samples,]
    taxcom_sample_ks = taxcom_sample_ks[! taxcom_sample_ks$sample %in% mhd_samples,]
    metrics_adval = metrics_adval[! metrics_adval$sample %in% mhd_samples,]
}else {
    focus="all" 
    print("Focus = all")
}

#####################################
# define groups of samples
healthy = rep("control", times = length(metadata[metadata$status == "ND", "sample"]))
names(healthy) = metadata[metadata$status == "ND", "sample"]

t2d = rep('T2D' , times = length(metadata[metadata$status == "T2D", "sample"]))
names(t2d) = metadata[metadata$status == "T2D", "sample"]

t1d = rep("T1D", times = length(metadata[metadata$status == "T1D", "sample"]))
names(t1d) = metadata[metadata$status == "T1D", "sample"]

diab = rep("diabetes", times = length(metadata[metadata$status == "T1D" | metadata$status == "T2D", "sample"]))
names(diab) = metadata[metadata$status == "T1D" | metadata$status == "T2D", "sample"]

t2d_tt = rep("T2D_metformin", times = length(metadata[metadata$status == "T2D" & metadata$type == "metformin+", "sample"]))
names(t2d_tt) = metadata[metadata$status == "T2D" & metadata$type == "metformin+", "sample"]

t2d_nott = rep("T2D_ctrl", times = length(metadata[metadata$status == "T2D" & metadata$type == "metformin-", "sample"]))
names(t2d_nott) = metadata[metadata$status == "T2D" & metadata$type == "metformin-", "sample"]

status_l1 = c(healthy, diab)
status_l1 = as.factor(status_l1)
status_l2 = c(healthy, t2d, t1d)
status_l2 = as.factor(status_l2)
status_l3 = c(healthy, t2d_tt, t2d_nott)
status_l3 = as.factor(status_l3)

#### colors 
statuses = c("control", "diabetes", "T1D", "T2D", "T2D_ctrl", "T2D_metformin")
# coles = viridis(length(statuses))
coles = c("#70B77E", "#00A5CF", "#453F3C", "#F75C03", "#F1C40F", "#7768AE") #, "#1C3144"
names(coles) = statuses
col_by_st = coles[as.character(status_l1)]
names(col_by_st) = names(status_l1)

taxcom_sample_ks = subset(taxcom_sample_ks, taxcom_sample_ks$expe == "addedvalue_com")

########################## Data treatment
# update factors levels in samples
taxcom_sample_ks$sample <- droplevels(taxcom_sample_ks$sample)
metrics_adval$sample <- droplevels(metrics_adval$sample)

metrics_adval$ks_ratio = metrics_adval$ks / metrics_adval$community_size
metrics_adval$es_ratio = metrics_adval$es / metrics_adval$community_size
metrics_adval$as_ratio = metrics_adval$as / metrics_adval$community_size
# create a df with all metrics and statuses
metrics = metrics_adval
names(metrics) = c("sample", "community_size", "ks", "es", "as", "ks_ratio", "es_ratio", "as_ratio")
metrics$status_l1 = status_l1[as.character(metrics$sample)]
metrics$status_l2 = status_l2[as.character(metrics$sample)]
metrics$status_l3 = status_l3[as.character(metrics$sample)]

# add category (init) and expe to taxcom_sample_init so that we can merge it with taxcom_sample_ks later
taxcom_sample_init$category = "init"
# make 2 tables, each with one expe value, merge them afterwards
t1 = taxcom_sample_init
t1$expe = "addedvalue_com"
# merge to taxcom_sample_ks
taxcom_sample_ks = rbind(taxcom_sample_ks, t1)

# add the initial community size
taxcom_sample_ks$total_init <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "community_size"]

# add total number of KS, ES, AS into the df
taxcom_sample_ks$total_KS <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "ks"]

taxcom_sample_ks$total_ES <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "es"]

taxcom_sample_ks$total_AS <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "as"]

# add the statuses
taxcom_sample_ks$status_l1 <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "status_l1"]
taxcom_sample_ks$status_l2 <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "status_l2"]
taxcom_sample_ks$status_l3 <- metrics[match(taxcom_sample_ks$sample, metrics$sample), "status_l3"]

# transform the occurrence into ratio depending on the considered category
taxcom_sample_ks$occurrence_ratio <- ifelse(
    taxcom_sample_ks$category == "AS", #test
    taxcom_sample_ks$occurrence/taxcom_sample_ks$total_AS,  #then
    ifelse( # else other test
        taxcom_sample_ks$category == "ES", #test
        taxcom_sample_ks$occurrence/taxcom_sample_ks$total_ES, # then
        ifelse(
            taxcom_sample_ks$category == "KS", #test
            taxcom_sample_ks$occurrence/taxcom_sample_ks$total_KS, # then
            ifelse(
                taxcom_sample_ks$category == "init", #test
                taxcom_sample_ks$occurrence/taxcom_sample_ks$total_init, # then
                NA
            )
        )
    )
)

taxtomerge = c("Cyanobacteria", "Verrucomicrobiota", "Spirochaetota", "Euryarchaeota", "Desulfobacterota_A", "Thermoplasmatota", "Synergistota")
taxcom_sample_ks$phylum_simp = ifelse(
            startsWith(as.character(taxcom_sample_ks$phylum), "Firmicutes"),
            "Firmicutes", 
            ifelse(as.character(taxcom_sample_ks$phylum) %in% taxtomerge,
                "Other",
                as.character(taxcom_sample_ks$phylum)
            )
) 
taxcom_sample_ks$phylum_simp = as.factor(taxcom_sample_ks$phylum_simp)
kept_tax = c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Other")
taxcom_sample_ks$phylum_simp = factor(taxcom_sample_ks$phylum_simp, levels = kept_tax)


taxcom_sample_ks$category = factor(taxcom_sample_ks$category, levels = c("init", "KS", "AS", "ES"))


# colors based on datasets
colds = c("#70B77E", "#F1C40F")
names(colds) = c("MHD", "SWE")
col_by_ds = colds[metadata$subset]
names(col_by_ds) = metadata$sample

#################### Boxplots
taxcolors = brewer.pal(n = 5, name = 'Dark2')
names(taxcolors) = kept_tax

for (status in c("status_l1", "status_l2", "status_l3")){
    current_df = taxcom_sample_ks#[taxcom_sample_ks$expe == "addedvalue_com",]
    p = ggplot(data=subset(current_df, !is.na(current_df[[status]])), aes(x = category, y = occurrence_ratio, fill = phylum_simp), show.legend = FALSE) +
    geom_boxplot() + 
    theme_classic() +
    facet_grid(~ .data[[status]], labeller = label_wrap_gen(width=10)) +
    # ggtitle("Taxonomic composition of keystone species") +
    scale_fill_manual(values=taxcolors) +
    # scale_y_continuous(trans='sqrt') +
    labs(y= "compositional ratio", x = "category", fill = "Phylum") +
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA),
        axis.text=element_text(size=15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title=element_text(size=15),
        strip.text.x=element_text(size=15)
        )
    pdf(paste0(res_dir, status, "_", "tax_compo.pdf"), width = 11, height = 7)
    print(p)
    dev.off()
}

#################### Tests
# # test for the effect of the category within a given category
# aov_res_status_effect <- data.frame(matrix(ncol = 6, nrow = 0))
# all_tests_status_effect <- list()
# colnames(aov_res_status_effect) <- c("expe", "phylum", "category", "pvalue_l1", "pvalue_l2", "pvalue_l3")
# for (expe in unique(as.character(taxcom_sample_ks$expe))) {
#     for (phylum in unique(as.character(taxcom_sample_ks$phylum_simp))) {
#         for (category in unique(as.character(taxcom_sample_ks$category))) {
#             datait = taxcom_sample_ks[taxcom_sample_ks$expe ==expe & taxcom_sample_ks$phylum_simp==phylum & taxcom_sample_ks$category==category, ]
#             fit1 = aov(occurrence_ratio ~ status_l1, datait)
#             fit2 = aov(occurrence_ratio ~ status_l2, datait)
#             fit3 = aov(occurrence_ratio ~ status_l3, datait)
#             all_tests_status_effect[[paste(expe, "_", phylum, "_", category, "_statusL1", sep='')]] <- fit1
#             all_tests_status_effect[[paste(expe, "_", phylum, "_", category, "_statusL2", sep='')]] <- fit2
#             all_tests_status_effect[[paste(expe, "_", phylum, "_", category, "_statusL3", sep='')]] <- fit3
#             aov_res_status_effect[nrow(aov_res_status_effect) + 1,] = c(expe, phylum, category, round(summary(fit1)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit2)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit3)[[1]][["Pr(>F)"]][1], digits = 3))        
#         }
#     }
# }


# # test for the effect of the category within a given status
# aov_res_category_effect <- data.frame(matrix(ncol = 5, nrow = 0))
# all_tests_category_effect <- list()
# colnames(aov_res_category_effect) <- c("expe", "phylum", "status_level", "status", "pvalue")
# for (expe in unique(as.character(taxcom_sample_ks$expe))) {
#     for (phylum in unique(as.character(taxcom_sample_ks$phylum_simp))) {
#         for (status_level in c("status_l1", "status_l2", "status_l3")) {
#             for (status in unique(na.omit(as.character(taxcom_sample_ks[[status_level]])))){
#                 datait = taxcom_sample_ks[taxcom_sample_ks$expe ==expe & taxcom_sample_ks$phylum_simp==phylum & taxcom_sample_ks[[status_level]]==status, ]
#                 datait = subset(datait, !is.na(datait[[status_level]]))
#                 fit1 = aov(occurrence_ratio ~ category, datait)
#                 all_tests_category_effect[[paste(expe, "_", phylum, "_", status_level, "_", status, sep='')]] <- fit1
#                 aov_res_category_effect[nrow(aov_res_category_effect) + 1,] = c(expe, phylum, status_level, status, round(summary(fit1)[[1]][["Pr(>F)"]][1], digits = 3))        
#             }
#         }
#     }
# }
