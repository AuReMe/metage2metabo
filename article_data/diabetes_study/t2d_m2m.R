# load packages
library(reshape2)
library(dplyr)
library(ape)
library(vegan)
library(viridis)
library(Rtsne)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

######### CONFIG ###########
analysis = "comscope" # can be "addedvalue" to study the cooperation potential or "comscope" for the full community scope
focus = "mhd" # can be "mhd" or "swe" or "all"
############################

if (analysis=="addedvalue"){
    # input files
    cscope_compo_file = 'addedvalue_long.tsv'
    cscope_onto_file = 'addedvalue_by_family.tsv'
} else{ 
    # input files
    cscope_compo_file = 'community_scopes_long.tsv'
    cscope_onto_file = 'community_scopes_by_family.tsv'
}

metadata_file = 'metadata.csv'
tax_file = 'mgs_specI.tax'
comm_size_file ="initial_community_size.tsv"
iscope_file="individual_scopes.tsv"
metrics_file = 'community_reduction_metrics.tsv'

res_dir = paste0('figures_results_', focus, "_", analysis, "/")
dir.create(file.path(res_dir), showWarnings = FALSE)


# we want to drop samples that had less than 1000 reads passing quality checks
samples_to_be_dropped = c("NG5425286", "NG5425485", "MH0035", "MH0037", "MH0038", "MH0047", "MH0054", "MH0057", "MH0058", "MH0064", "MH0065", "MH0067", "MH0073", "MH0077", "MH0082", "MH0084")

# read input files
comm_size = as.data.frame(read.table(comm_size_file, header=TRUE, sep="\t"))
names(comm_size) = c("sample", "comsize")
comm_size$sample = as.character(comm_size$sample)
cscope_compo <- as.data.frame(read.table(cscope_compo_file, header=FALSE, sep="\t"))
colnames(cscope_compo) = c("compound", "sample", "occurrence")
# there might be duplicated rows (compounds producible in multiple compartments originally)
cscope_compo = distinct(cscope_compo, .keep_all = TRUE)

cscope_onto <- as.data.frame(read.table(cscope_onto_file, header=FALSE, sep="\t"))
colnames(cscope_onto) = c("sample", "family", "occurrence")

metrics = as.data.frame(read.table(metrics_file, header=TRUE, sep="\t"))

# create matrices
compo = acast(cscope_compo, sample~compound, value.var="occurrence", fill=0)
onto = acast(cscope_onto, sample~family, value.var="occurrence", fill=0)


# read metadata
metadata = as.data.frame(read.table(metadata_file, header=TRUE, sep=","))
tax = as.data.frame(read.table(tax_file, header=FALSE, sep="\t"))
names(tax) = c("genome","K","P","C","O","F","G","S")
# remove metadata samples that are missing
metadata = metadata[! metadata$sample %in% samples_to_be_dropped, ]

# read iscope file
iscope_long = as.data.frame(read.table(iscope_file, header=FALSE, sep="\t"))
names(iscope_long) = c("compound", "genome", "occurrence")
# there might be duplicated rows (compounds producible in multiple compartments originally)
iscope_long = distinct(iscope_long, .keep_all = TRUE)
iscope = acast(iscope_long, genome~compound, value.var="occurrence", fill=0)

mhd_samples = metadata[metadata$subset == 'MHD', 'sample']
swe_samples = metadata[metadata$subset == 'SWE', 'sample']

if (focus=="mhd"){
    # drop the SWE samples
    cscope_onto = cscope_onto[! cscope_onto$sample %in% swe_samples,]
    cscope_compo = cscope_compo[! cscope_compo$sample %in% swe_samples,]
    compo = compo[! rownames(compo) %in% swe_samples,]
    onto = onto[! rownames(onto) %in% swe_samples,]
    metadata = metadata[! metadata$sample %in% swe_samples,]
} else if (focus=="swe"){
    # drop the MHD samples
    cscope_onto = cscope_onto[! cscope_onto$sample %in% mhd_samples,]
    cscope_compo = cscope_compo[! cscope_compo$sample %in% mhd_samples,]
    compo = compo[! rownames(compo) %in% mhd_samples,]
    onto = onto[! rownames(onto) %in% mhd_samples,]
    metadata = metadata[! metadata$sample %in% mhd_samples,]
}else {
    focus="all" 
    print("Focus = all")
    write.table(x=compo, file=paste0(res_dir, "compo.tsv"), sep = '\t', quote = FALSE, col.names = NA)
    write.table(x=onto, file=paste0(res_dir, "onto.tsv"), sep = '\t', quote = FALSE, col.names = NA)
}

# in any case, drop the samples that had not enough reads
cscope_onto = cscope_onto[! cscope_onto$sample %in% samples_to_be_dropped,]
compo = compo[! rownames(compo) %in% samples_to_be_dropped,]
onto = onto[! rownames(onto) %in% samples_to_be_dropped,]
comm_size = comm_size[! comm_size$sample %in% samples_to_be_dropped,]
metrics = metrics[metrics$sample %in% rownames(onto),]

# remove columns of compo and onto that end up with only zeros since we removed samples
compo = compo[,-(which(colSums(compo)==0))]
onto = onto[,-(which(colSums(onto)==0))]


# retrieve tested values (compounds and families)
all_families = colnames(onto)
all_compounds = colnames(compo)

# ratio of KS out of total community_size
metrics$ks_ratio = metrics$ks / metrics$community_size
metrics$es_ratio = metrics$es / metrics$community_size
metrics$as_ratio = metrics$as / metrics$community_size

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

# add status to comsize df
comm_size$status_l1 = status_l1[as.character(comm_size$sample)]
comm_size$status_l2 = status_l2[as.character(comm_size$sample)]
comm_size$status_l3 = status_l3[as.character(comm_size$sample)]

# add the status to the long dataframes
cscope_compo$status_l1 = status_l1[as.character(cscope_compo$sample)]
cscope_compo$status_l2 = status_l2[as.character(cscope_compo$sample)]
cscope_compo$status_l3 = status_l3[as.character(cscope_compo$sample)]

cscope_onto$status_l1 = status_l1[as.character(cscope_onto$sample)]
cscope_onto$status_l2 = status_l2[as.character(cscope_onto$sample)]
cscope_onto$status_l3 = status_l3[as.character(cscope_onto$sample)]

# size of cscope by sample and status
cscope_size = as.data.frame(rowSums(onto))
names(cscope_size) = c("cscope_sz")
cscope_size$sample = rownames(cscope_size)
cscope_size$status_l1 = status_l1[as.character(cscope_size$sample)]
cscope_size$status_l2 = status_l2[as.character(cscope_size$sample)]
cscope_size$status_l3 = status_l3[as.character(cscope_size$sample)]
cscope_size$comsize = comm_size[comm_size$sample %in% cscope_size$sample, "comsize"]
cscope_size$cohort = metadata[match(cscope_size$sample, metadata$sample), "subset"]


onto_l1 = aggregate(cscope_onto$occurrence, by=list(status=cscope_onto$status_l1, family=cscope_onto$family), FUN=sum, na.rm=TRUE)
# turn into matrix
onto_l1 = acast(onto_l1, family~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(onto_l1)){
    onto_l1[,col] <- onto_l1[,col] / summary(status_l1)[col]
}
onto_l1_long = melt(onto_l1)
names(onto_l1_long) = c("family", "status", "mean_occurrence")

onto_l2 = aggregate(cscope_onto$occurrence, by=list(status=cscope_onto$status_l2, family=cscope_onto$family), FUN=sum, na.rm=TRUE)
# turn into matrix
onto_l2 = acast(onto_l2, family~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(onto_l2)){
    onto_l2[,col] <- onto_l2[,col] / summary(status_l2)[col]
}
onto_l2_long = melt(onto_l2)
names(onto_l2_long) = c("family", "status", "mean_occurrence")

onto_l3 = aggregate(cscope_onto$occurrence, by=list(status=cscope_onto$status_l3, family=cscope_onto$family), FUN=sum, na.rm=TRUE)
# turn into matrix
onto_l3 = acast(onto_l3, family~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(onto_l3)){
    onto_l3[,col] <- onto_l3[,col] / summary(status_l3)[col]
}
onto_l3_long = melt(onto_l3)
names(onto_l3_long) = c("family", "status", "mean_occurrence")

compo_l1 = aggregate(cscope_compo$occurrence, by=list(status=cscope_compo$status_l1, compound=cscope_compo$compound), FUN=sum, na.rm=TRUE)
# turn into matrix
compo_l1 = acast(compo_l1, compound~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(compo_l1)){
    compo_l1[,col] <- compo_l1[,col] / summary(status_l1)[col]
}
compo_l1_long = melt(compo_l1)
names(compo_l1_long) = c("compound", "status", "mean_occurrence")


compo_l2 = aggregate(cscope_compo$occurrence, by=list(status=cscope_compo$status_l2, compound=cscope_compo$compound), FUN=sum, na.rm=TRUE)
# turn into matrix
compo_l2 = acast(compo_l2, compound~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(compo_l2)){
    compo_l2[,col] <- compo_l2[,col] / summary(status_l2)[col]
}
compo_l2_long = melt(compo_l2)
names(compo_l2_long) = c("compound", "status", "mean_occurrence")

compo_l3 = aggregate(cscope_compo$occurrence, by=list(status=cscope_compo$status_l3, compound=cscope_compo$compound), FUN=sum, na.rm=TRUE)
# turn into matrix
compo_l3 = acast(compo_l3, compound~status, value.var="x", fill=0)
# divide by number of samples in each status category
for (col in colnames(compo_l3)){
    compo_l3[,col] <- compo_l3[,col] / summary(status_l3)[col]
}
compo_l3_long = melt(compo_l3)
names(compo_l3_long) = c("compound", "status", "mean_occurrence")

############ TESTS ##############
# AOV tests init comm size  and commscope size
init_size_aov_l1 = aov(comsize ~ status_l1, cscope_size)
init_size_aov_l2 = aov(comsize ~ status_l2, cscope_size)
init_size_aov_l3 = aov(comsize ~ status_l3, cscope_size)

summary(init_size_aov_l1)
TukeyHSD(init_size_aov_l1)
summary(init_size_aov_l2)
TukeyHSD(init_size_aov_l2)
summary(init_size_aov_l3)
TukeyHSD(init_size_aov_l3)

group_by(cscope_size, status_l2) %>%
  summarise(
    count = n(),
    mean = mean(comsize, na.rm = TRUE),
    sd = sd(comsize, na.rm = TRUE),
  )

cscope_size_aov_l1 = aov(cscope_sz ~ status_l1, cscope_size)
cscope_size_aov_l2 = aov(cscope_sz ~ status_l2, cscope_size)
cscope_size_aov_l3 = aov(cscope_sz ~ status_l3, cscope_size)

summary(cscope_size_aov_l1)
TukeyHSD(cscope_size_aov_l1)
summary(cscope_size_aov_l2)
TukeyHSD(cscope_size_aov_l2)
summary(cscope_size_aov_l3)
TukeyHSD(cscope_size_aov_l3)

group_by(cscope_size, status_l2) %>%
  summarise(
    count = n(),
    mean = mean(cscope_sz, na.rm = TRUE),
    sd = sd(cscope_sz, na.rm = TRUE),
  )


# AOV TEST compounds
aov_per_family <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(aov_per_family) <- c("family", "pvalue_l1", "pvalue_l2", "pvalue_l3")#, "pvalue_cohort")
for (fam in c(as.character(all_families))){
    datait = select(as.data.frame(onto), fam)
    datait$sample = rownames(datait)
    datait$status_l1 = status_l1[as.character(datait$sample)]
    datait$status_l2 = status_l2[as.character(datait$sample)]
    datait$status_l3 = status_l3[as.character(datait$sample)]
    datait$cohort = metadata[match(datait$sample,metadata$sample), "subset"]
    fit1 = aov(datait[[fam]] ~ status_l1, datait)
    fit2 = aov(datait[[fam]] ~ status_l2, datait)
    fit3 = aov(datait[[fam]] ~ status_l3, datait)
    # fitc = aov(datait[[fam]] ~ cohort, datait)
    aov_per_family[nrow(aov_per_family) + 1,] = c(fam, round(summary(fit1)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit2)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit3)[[1]][["Pr(>F)"]][1], digits = 3))#, round(summary(fitc)[[1]][["Pr(>F)"]][1], digits = 3))
}

aov_per_compound <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(aov_per_compound) <- c("compound", "pvalue_l1", "pvalue_l2", "pvalue_l3")
for (fam in c(as.character(all_compounds))){
    datait = select(as.data.frame(compo), fam)
    datait$sample = rownames(datait)
    datait$status_l1 = status_l1[as.character(datait$sample)]
    datait$status_l2 = status_l2[as.character(datait$sample)]
    datait$status_l3 = status_l3[as.character(datait$sample)]
    fit1 = aov(datait[[fam]] ~ status_l1, datait)
    fit2 = aov(datait[[fam]] ~ status_l2, datait)
    fit3 = aov(datait[[fam]] ~ status_l3, datait)
    aov_per_compound[nrow(aov_per_compound) + 1,] = c(fam, round(summary(fit1)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit2)[[1]][["Pr(>F)"]][1], digits = 3), round(summary(fit3)[[1]][["Pr(>F)"]][1], digits = 3))
}

# aov_per_compound[grep("utyr", aov_per_compound$compound, ignore.case=TRUE), ]
# aov_per_family[grep("lip", aov_per_family$family, ignore.case=TRUE), ]
# compo_l3[grep("utyr", rownames(compo_l3), ignore.case = TRUE), ]
# compo_l2[grep("utyr", rownames(compo_l2), ignore.case = TRUE), ]
# compo_l1[grep("utyr", rownames(compo_l1), ignore.case = TRUE), ]
# onto_l3[grep("lip", rownames(onto_l3), ignore.case = TRUE), ]
# onto_l2[grep("lip", rownames(onto_l2), ignore.case = TRUE), ]
# onto_l1[grep("lip", rownames(onto_l1), ignore.case = TRUE), ]

#### colors 
statuses = c("control", "diabetes", "T1D", "T2D", "T2D_ctrl", "T2D_metformin")
# coles = viridis(length(statuses))
coles = c("#70B77E", "#00A5CF", "#453F3C", "#F75C03", "#F1C40F", "#7768AE")#, "#1C3144"
names(coles) = statuses
col_by_st = coles[as.character(status_l1)]
names(col_by_st) = names(status_l1)

# colors based on datasets
colds = c("#70B77E", "#F1C40F")
names(colds) = c("MHD", "SWE")
col_by_ds = colds[metadata$subset]
names(col_by_ds) = metadata$sample

############ figures
for (status in c("status_l1", "status_l2", "status_l3")){
    current_colors = coles[as.character(unique(cscope_size[[status]]))]
    p <- ggplot(data = subset(cscope_size, !is.na(cscope_size[[status]])), aes(x=.data[[status]], y=cscope_sz, fill=.data[[status]])) + geom_violin() + theme_classic() + scale_fill_manual(values=current_colors) + geom_boxplot(width=0.1) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=18), axis.title=element_text(size=21),) + labs(y= "community scope size", x = "disease status")

    q <- ggplot(data = subset(cscope_size, !is.na(cscope_size[[status]])), aes(x=.data[[status]], y=comsize, fill=.data[[status]])) + geom_violin() + theme_classic() + scale_fill_manual(values=current_colors) + geom_boxplot(width=0.1) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=18), axis.title=element_text(size=21),) + labs(y= "initial community size", x = "disease status")

    pdf(paste0(res_dir, status, "_", "violin_plot_comscope_size.pdf"))
    print(p)
    dev.off()
    pdf(paste0(res_dir, status, "_", "violin_plot_comsize.pdf")) 
    print(q)
    dev.off()
}


cohortcols = c("#7f5da2", "#f4c84f")
names(cohortcols) = c("MHD", "SWE")
if (focus == "all"){
    if (analysis == "comscope"){
        ylabtitle = "community scope size"
    }else if (analysis == "addedvalue"){
        ylabtitle = "cooperation potential size"
    }
    r <- ggplot(data = cscope_size, aes(x=.data[["cohort"]], y=comsize, fill=.data[["cohort"]])) + geom_violin() + theme_classic() + scale_fill_manual(values=cohortcols) + geom_boxplot(width=0.1) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=15),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12),
    axis.title=element_text(size=15),
    strip.text.x=element_text(size=12)) + labs(y= "community size", x = "cohort")
    pdf(paste0(res_dir, "cohort", "_", "violin_plot_comsize.pdf"))
    print(r)
    dev.off()
    s <- ggplot(data = cscope_size, aes(x=.data[["cohort"]], y=cscope_sz, fill=.data[["cohort"]])) + geom_violin() + theme_classic() + scale_fill_manual(values=cohortcols) + geom_boxplot(width=0.1) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=15),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12),
    axis.title=element_text(size=15),
    strip.text.x=element_text(size=12)) + labs(y= ylabtitle, x = "cohort")
    pdf(paste0(res_dir, "cohort", "_", "violin_plot_scope_size.pdf"))
    print(s)
    dev.off()
}

statuses = c("control", "diabetes", "T1D", "T2D", "T2D_ctrl", "T2D_metformin")

################## correlation analysis (spearman)
kept3 = c("Acids", "Alcohols", "Aldehydes-Or-Ketones", "All-Amines", "All-Amino-Acids", "All-Carbohydrates", "All-Nucleosides", "Amides", "Antibiotics", "Aromatics", "Cofactors", "Esters", "Hormones", "Lipids", "Nitrogen-Molecular-Entities", "Organic-heterocyclic-compound", "ORGANOSULFUR", "Secondary-Metabolites", "Others")

onto_red = onto[,colnames(onto) %in% kept3]
onto_redf = as.data.frame(onto_red)
onto_redf$status = status_l1[rownames(onto_redf)]
# onto_red = t(onto_red)
agg_onto = aggregate(onto_redf[,colnames(onto_redf) %in% c(as.character(kept3))], list(onto_redf$status), mean, na.rm=TRUE)
rownames(agg_onto) <- agg_onto[,1]
agg_onto <- agg_onto[,-1]

cormat <- round(cor(onto_red, method="spearman"),2)
# head(cormat)
# head(melted_cormat)
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

pdf(paste0(res_dir, "corr_spearman_ontofam.pdf")) 
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
geom_tile(color = "white")+
scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
midpoint = 0, limit = c(-1,1), space = "Lab",
name="Spearman\nCorrelation") +
theme_minimal()+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, 
size = 12, hjust = 1))+
coord_fixed()
dev.off()

################### boxplot of families by subset
kept4 = c("Acids", "All-Amino-Acids", "All-Carbohydrates", "All-Nucleosides", "Aromatics", "Lipids", "Secondary-Metabolites")

kept5 = c("Acids", "Aldehydes-Or-Ketones", "All-Carbohydrates",  "Aromatics", "Lipids", "Nitrogen-Molecular-Entities", "Organic-heterocyclic-compound", "ORGANOSULFUR", "Secondary-Metabolites", "Others")

if (analysis =="comscope"){
    kept = kept4
    cscope_onto$cohort = metadata[match(cscope_onto$sample, metadata$sample), "subset"]
    cohort_cscope = cscope_onto[cscope_onto$family %in% kept,]
    cohort_cscope$family = gsub("-", " ", cohort_cscope$family)
    p = ggplot(data=cohort_cscope, aes(x = cohort, y = occurrence, fill = cohort), show.legend = FALSE) +
            geom_boxplot() + 
            theme_classic() +
            facet_grid(~ .data[["family"]], labeller = label_wrap_gen(width=10)) +
            #ggtitle("Community scope composition by cohort") +
            scale_fill_manual(values=cohortcols) +
            # scale_y_continuous(trans='sqrt') +
            labs(y= "occurrence of family in community scopes", x = "metabolites family", fill = "Metabolites family") +
            theme(panel.grid = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none",
                panel.border = element_rect(colour = "black", fill=NA),
                axis.text=element_text(size=15),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12),
                axis.title=element_text(size=15),
                strip.text.x=element_text(size=15)
                )
    pdf(paste0(res_dir, analysis, "_", focus, "_cohort_", "onto_some_families.pdf"), width = 18, height = 7)
    print(p)
    dev.off()
}else if (analysis == "addedvalue"){
    kept = kept5
    cscope_onto$cohort = metadata[match(cscope_onto$sample, metadata$sample), "subset"]
    cohort_cscope = cscope_onto[cscope_onto$family %in% kept,]
    cohort_cscope$family = gsub("-", " ", cohort_cscope$family)
    cohort_cscope$family = gsub("ORGANOSULFUR", "Organosulfur", cohort_cscope$family)
    p = ggplot(data=cohort_cscope, aes(x = cohort, y = occurrence, fill = cohort), show.legend = FALSE) +
            geom_boxplot() + 
            theme_classic() +
            facet_grid(~ .data[["family"]], labeller = label_wrap_gen(width=10)) +
            #ggtitle("Community scope composition by cohort") +
            scale_fill_manual(values=cohortcols) +
            # scale_y_continuous(trans='sqrt') +
            labs(y= "occurrence of family in added value of cooperation", x = "metabolites family", fill = "Metabolites family") +
            theme(panel.grid = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none",
                panel.border = element_rect(colour = "black", fill=NA),
                axis.text=element_text(size=15),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12),
                axis.title=element_text(size=15),
                strip.text.x=element_text(size=12)
                )
    pdf(paste0(res_dir, analysis, "_", focus, "_cohort_", "onto_some_families.pdf"), width = 18, height = 7)
    print(p)
    dev.off()
}

################## Heatmap of compo and onto matrices
pdf(paste0(res_dir, "heatmap_compo_coloredby_dataset.pdf")) 
heatmap(compo, RowSideColors = col_by_ds[rownames(compo)], scale = "none")
legend("topright", legend=names(colds), pch=16, col=colds)
dev.off()
pdf(paste0(res_dir, "heatmap_onto_coloredby_dataset.pdf")) 
heatmap(onto_red, RowSideColors = col_by_ds[rownames(onto_red)])
legend("topright", legend=names(colds), pch=16, col=colds)
dev.off()

boolcolors = c("white", "orange")
names(boolcolors) = c(0,1)
library(gplots)
for (status in c("status_l1", "status_l2", "status_l3")){
    if (status == "status_l1"){st <- status_l1}
    else if (status == "status_l2") {st <- status_l2}
    else {st <- status_l3}
    col_by_st = coles[as.character(st)]
    names(col_by_st) = names(st)
    current_colors = coles[as.character(unique(st))]
    pdf(paste0(res_dir, status, "_heatmap_compo.pdf")) 
    heatmap(compo, RowSideColors = col_by_st[rownames(compo)], scale = "none", cexRow=0.4, cexCol=0.6, col = boolcolors)
    legend("topright", legend=names(current_colors), pch=16, col=current_colors)
    dev.off()
    pdf(paste0(res_dir, status, "_heatmap_onto.pdf")) 
    heatmap.2(onto_red, RowSideColors = col_by_st[rownames(onto_red)], cexRow=0.4, cexCol=0.6, srtCol = 45, adjCol = c(1,0), offsetCol=1, col    = brewer.pal(9, "OrRd"),density.info="none", trace="none")
    legend("topright", legend=names(current_colors), pch=16, col=current_colors)
    dev.off()
}
