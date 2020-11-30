# load packages
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

######### CONFIG ###########
focus = "mhd" # choose between "mhd" "swe" and "all"
############################

metadata_file = 'metadata.csv'
metrics_file_av = 'community_reduction_metrics.tsv'

res_dir = paste0('figures_results_', focus, "_sizes", "/")
dir.create(file.path(res_dir), showWarnings = FALSE)

# we want to drop samples that had less than 1000 reads passing quality checks
samples_to_be_dropped = c("NG5425286", "NG5425485", "MH0035", "MH0037", "MH0038", "MH0047", "MH0054", "MH0057", "MH0058", "MH0064", "MH0065", "MH0067", "MH0073", "MH0077", "MH0082", "MH0084")

# read input files
metrics_adval = as.data.frame(read.table(metrics_file_av, header=TRUE, sep="\t"))
metadata = as.data.frame(read.table(metadata_file, header=TRUE, sep=","))

# in any case, drop the samples that had not enough reads
metrics_adval = metrics_adval[! metrics_adval$sample %in% samples_to_be_dropped,]
metadata = metadata[! metadata$sample %in% samples_to_be_dropped, ]

mhd_samples = metadata[metadata$subset == 'MHD', 'sample']
swe_samples = metadata[metadata$subset == 'SWE', 'sample']

if (focus=="mhd"){
    # drop the SWE samples
    metrics_adval = metrics_adval[! metrics_adval$sample %in% swe_samples,]
    metadata = metadata[! metadata$sample %in% swe_samples,]
} else if (focus=="swe"){
    # drop the MHD samples
    metrics_adval = metrics_adval[! metrics_adval$sample %in% mhd_samples,]
    metadata = metadata[! metadata$sample %in% mhd_samples,]
}else {
    focus="all" 
    print("Focus = all")
}

# ratio of KS out of total community_size
metrics_adval$ks_ratio = metrics_adval$ks / metrics_adval$community_size
metrics_adval$es_ratio = metrics_adval$es / metrics_adval$community_size
metrics_adval$as_ratio = metrics_adval$as / metrics_adval$community_size

# define groups of samples
healthy = rep("control", times = length(metadata[metadata$status == "ND", "sample"]))
names(healthy) = metadata[metadata$status == "ND", "sample"]
# names(names(healthy)) = c(rep("control", length(healthy)))

t2d = rep('T2D' , times = length(metadata[metadata$status == "T2D", "sample"]))
names(t2d) = metadata[metadata$status == "T2D", "sample"]
# names(t2d) = c(rep("T2D", length(t2d)))

t1d = rep("T1D", times = length(metadata[metadata$status == "T1D", "sample"]))
names(t1d) = metadata[metadata$status == "T1D", "sample"]
# names(t1d) = c(rep("T1D", length(t1d)))

diab = rep("diabetes", times = length(metadata[metadata$status == "T1D" | metadata$status == "T2D", "sample"]))
names(diab) = metadata[metadata$status == "T1D" | metadata$status == "T2D", "sample"]
# names(diab) = c(rep("diabetes", length(diab)))

t2d_tt = rep("T2D_metformin", times = length(metadata[metadata$status == "T2D" & metadata$type == "metformin+", "sample"]))
names(t2d_tt) = metadata[metadata$status == "T2D" & metadata$type == "metformin+", "sample"]
# names(t2d_tt) = c(rep("T2D_metformin", length(t2d_tt)))

t2d_nott = rep("T2D_ctrl", times = length(metadata[metadata$status == "T2D" & metadata$type == "metformin-", "sample"]))
names(t2d_nott) = metadata[metadata$status == "T2D" & metadata$type == "metformin-", "sample"]
# names(t2d_nott) = c(rep("T2D_ctrl", length(t2d_nott)))

status_l1 = c(healthy, diab)
status_l1 = as.factor(status_l1)
status_l2 = c(healthy, t2d, t1d)
status_l2 = as.factor(status_l2)
status_l3 = c(healthy, t2d_tt, t2d_nott)
status_l3 = as.factor(status_l3)

# create a df with all metrics and statuses
metrics = metrics_adval
names(metrics) = c("sample", "community_size", "ks", "es", "as", "ks_ratio", "es_ratio", "as_ratio")
metrics$status_l1 = status_l1[as.character(metrics$sample)]
metrics$status_l2 = status_l2[as.character(metrics$sample)]
metrics$status_l3 = status_l3[as.character(metrics$sample)]


#### colors 
statuses = c("control", "diabetes", "T1D", "T2D", "T2D_ctrl", "T2D_metformin")
# coles = viridis(length(statuses))
coles = c("#70B77E", "#00A5CF", "#453F3C", "#F75C03", "#F1C40F", "#7768AE", "#1C3144")
names(coles) = statuses
col_by_st = coles[as.character(status_l1)]
names(col_by_st) = names(status_l1)

# colors based on datasets
colds = c("#70B77E", "#F1C40F")
names(colds) = c("MHD", "SWE")
col_by_ds = colds[metadata$subset]
names(col_by_ds) = metadata$sample

##################################################################################
# Stats on KS, ES, AS sizes for addedvalue targets and full targets (addedvalue + community scope)
# ksafit1 = aov(ks ~ status_l1, metrics)
# ksafit2 = aov(ks ~ status_l2, metrics)
# ksafit3 = aov(ks ~ status_l3, metrics)
# esafit1 = aov(es ~ status_l1, metrics)
# esafit2 = aov(es ~ status_l2, metrics)
# esafit3 = aov(es ~ status_l3, metrics)
# asafit1 = aov(as ~ status_l1, metrics)
# asafit2 = aov(as ~ status_l2, metrics)
# asafit3 = aov(as ~ status_l3, metrics)

# summary(ksafit1) # *
# summary(ksafit2)
# summary(ksafit3)
# summary(esafit1)
# summary(esafit2)
# summary(esafit3)
# summary(asafit1) # *
# summary(asafit2)
# summary(asafit3)
# # if yes, Tukey multiple pairwise-comparisons
# TukeyHSD(ksafit1) # *
# TukeyHSD(ksafit2)
# TukeyHSD(ksafit3)
# TukeyHSD(esafit1)
# TukeyHSD(esafit2)
# TukeyHSD(esafit3)
# TukeyHSD(asafit1) # *
# TukeyHSD(asafit2)
# TukeyHSD(asafit3)


# Plots
plots_list = list()
for (status in c("status_l1", "status_l2", "status_l3")){
    current_colors = coles[as.character(unique(metrics[[status]]))]
    current_colors<-current_colors[!is.na(current_colors)]
    for (var in c("community_size", "ks", "ks_ratio", "es", "es_ratio", "as", "as_ratio")){
        current_plot = ggplot(subset(metrics, !is.na(metrics[[status]])), aes(x=.data[[status]], y=.data[[var]], fill=.data[[status]])) + geom_violin() + theme_classic() + scale_fill_manual(values=current_colors) + geom_boxplot(width=0.1) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), axis.text=element_text(size=7),
        axis.title=element_text(size=9), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA)) + labs(x = "disease status")
        plots_list[[paste("plot_", status, "_", var, sep='')]] <- current_plot
    }
}

all_plots = ggarrange(plotlist=plots_list, nrow = 3, ncol = 3)

pdf(paste0(res_dir, "violin_plots_commsizes_and_ratios_", focus, ".pdf"))
all_plots
dev.off()
