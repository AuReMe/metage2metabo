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
focus = "all" # can be "mhd" or "swe" or "all"
############################

metadata_file = 'metadata.csv'
tax_file = 'mgs_specI.tax'
abmatrix_file ="abundance.mat"
cscope_compo_file = 'community_scopes_long.tsv'
iscope_file="individual_scopes.tsv"

res_dir = paste0('figures_results_diversity_', focus, "/")
dir.create(file.path(res_dir), showWarnings = FALSE)

# we want to drop samples that had less than 1000 reads passing quality checks
samples_to_be_dropped = c("NG5425286", "NG5425485", "MH0035", "MH0037", "MH0038", "MH0047", "MH0054", "MH0057", "MH0058", "MH0064", "MH0065", "MH0067", "MH0073", "MH0077", "MH0082", "MH0084")

# read metadata
metadata = as.data.frame(read.table(metadata_file, header=TRUE, sep=","))
tax = as.data.frame(read.table(tax_file, header=FALSE, sep="\t"))
names(tax) = c("genome","K","P","C","O","F","G","S")
# remove metadata samples that are missing
metadata = metadata[! metadata$sample %in% samples_to_be_dropped, ]
# gets statuses
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

# read abundance matrix
mx = as.matrix(read.table(abmatrix_file, header=TRUE, row.names=1, sep="\t"))
print(paste0("Data loaded: ", nrow(mx), " bacteria, ", ncol(mx), " samples. "))
# remove unknown bact
mx = mx[!rownames(mx) %in% c("?", "?;?;?;?;?;?","?;?;?;?;?","?;?;?;?","?;?;?","?;?"),]
# print how many removed bact
print(paste0("After removing unknown bacteria: ", nrow(mx) , " bacteria, ", ncol(mx), " samples. "))
# remove bact that were never observed
mx = mx[rowSums(mx > 0) > 0,]
# print how many removed bact
print(paste0("After removing non-observed bacteria: ", nrow(mx), " bacteria, ", ncol(mx), " samples. "))
# remove samples that were not observed
mx = mx[,colSums(mx > 0) > 0]
# print how many removed samples
print(paste0("After removing unobserved samples: ", nrow(mx), " bacteria, ", ncol(mx), " samples. "))
# drop samples that had less than 1000 reads passing quality checks
mx = mx[,! colnames(mx) %in% samples_to_be_dropped]
# normalize by sum of all reads in a sample (i.e. column)
mx_norm = sweep(mx,2,colSums(mx),"/")

# MHD and SWE samples 
mhd_samples = metadata[metadata$subset == 'MHD', 'sample']
swe_samples = metadata[metadata$subset == 'SWE', 'sample']

# diversity
shannon_div = diversity(t(mx_norm))
# average diversity by group
tapply(diversity(t(mx_norm)), metadata$subset, mean)
# richness
richness = specnumber(t(mx_norm), MARGIN = 1)
richness_gp = specnumber(t(mx_norm), metadata$subset, MARGIN = 1)

df <- as.data.frame(shannon_div)
df = as.data.frame(df)
df <- tibble::rownames_to_column(df, "ID")
df$richness <- richness[df$ID]
df$subset <- metadata[match(df$ID, metadata$sample), "subset"]
df$status_l1 <- status_l1[df$ID]
df$status_l2 <- status_l2[df$ID]
df$status_l3 <- status_l3[df$ID]

#### colors 
statuses = c("control", "diabetes", "T1D", "T2D", "T2D_ctrl", "T2D_metformin")
# coles = viridis(length(statuses))
coles = c("#70B77E", "#00A5CF", "#453F3C", "#F75C03", "#F1C40F", "#7768AE") #, "#1C3144"
names(coles) = statuses
col_by_st = coles[as.character(status_l1)]
names(col_by_st) = names(status_l1)
# colors based on datasets
colds = c("#70B77E", "#F1C40F")
names(colds) = c("MHD", "SWE")
col_by_ds = colds[metadata$subset]
names(col_by_ds) = metadata$sample


# metabolic landscape of samples (with cooperation)
cscope_compo <- as.data.frame(read.table(cscope_compo_file, header=FALSE, sep="\t"))
colnames(cscope_compo) = c("compound", "sample", "occurrence")
# there might be duplicated rows (compounds producible in multiple compartments originally)
cscope_compo = distinct(cscope_compo, .keep_all = TRUE)
compo = acast(cscope_compo, sample~compound, value.var="occurrence", fill=0)
# same without cooperation (individual scopes)
iscope_long = as.data.frame(read.table(iscope_file, header=FALSE, sep="\t"))
names(iscope_long) = c("compound", "genome", "occurrence")
# there might be duplicated rows (compounds producible in multiple compartments originally)
iscope_long = distinct(iscope_long, .keep_all = TRUE)
iscope = acast(iscope_long, genome~compound, value.var="occurrence", fill=0)

if (focus=="mhd"){
    # drop the SWE samples
    compo = compo[! rownames(compo) %in% swe_samples,]
    df = df[! df$ID %in% swe_samples,]
} else if (focus=="swe"){
    # drop the MHD samples
    compo = compo[! rownames(compo) %in% mhd_samples,]
    df = df[! df$ID %in% mhd_samples,]
}else {
    focus="all" 
    print("Focus = all")
}

pdf(paste0(res_dir, focus, "_hist_shannon_col_by_ds.pdf"), width = 9, height = 7 ) 
    s = ggplot(df, aes(x = shannon_div)) +
        geom_histogram(aes(color = subset, fill = subset), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=colds) +
        scale_fill_manual(values=colds)
    print(s)
dev.off()
        
pdf(paste0(res_dir, focus, "_hist_richness_col_by_ds.pdf"), width = 9, height = 7 ) 
    r = ggplot(df, aes(x = richness)) +
        geom_histogram(aes(color = subset, fill = subset), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=colds) +
        scale_fill_manual(values=colds)
    print(r)
dev.off()

pdf(paste0(res_dir, focus, "_hist_shannon.pdf"), width = 9, height = 7 ) 
    s = ggplot(df, aes(x = shannon_div)) +
        geom_histogram(aes(), color="black", fill="white",
                        position = "identity", bins = 20, alpha = 0.4) +
        theme_classic() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=18), axis.title=element_text(size=21),) + labs(y= "count", x = "shannon diversity index")
    print(s)
dev.off()
        
pdf(paste0(res_dir, focus, "_hist_richness.pdf"), width = 9, height = 7 ) 
    r = ggplot(df, aes(x = richness)) +
        geom_histogram(aes(), color="black", fill="white",
                        position = "identity", bins = 20, alpha = 0.4) +
        theme_classic() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), axis.text=element_text(size=18), axis.title=element_text(size=21),) + labs(y= "count", x = "richness")
    print(r)
dev.off()

pdf(paste0(res_dir, focus, "_hist_richness_col_by_st1.pdf"), width = 9, height = 7 ) 
    r1 = ggplot(df, aes(x = richness)) +
        geom_histogram(aes(color = status_l1, fill = status_l1), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=coles) +
        scale_fill_manual(values=coles)
    print(r1)
dev.off()

pdf(paste0(res_dir, focus, "_hist_richness_col_by_st2.pdf"), width = 9, height = 7 ) 
    r2 = ggplot(df, aes(x = richness)) +
        geom_histogram(aes(color = status_l2, fill = status_l2), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=coles) +
        scale_fill_manual(values=coles)
    print(r2)
dev.off()

pdf(paste0(res_dir, focus, "_hist_shannon_col_by_st1.pdf"), width = 9, height = 7 ) 
    s1 = ggplot(df, aes(x = shannon_div)) +
        geom_histogram(aes(color = status_l1, fill = status_l1), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=coles) +
        scale_fill_manual(values=coles)
    print(s1)
dev.off()

pdf(paste0(res_dir, focus, "_hist_shannon_col_by_st2.pdf"), width = 9, height = 7 ) 
    s2 = ggplot(df, aes(x = shannon_div)) +
        geom_histogram(aes(color = status_l2, fill = status_l2), 
                        position = "identity", bins = 30, alpha = 0.4) +
        scale_color_manual(values=coles) +
        scale_fill_manual(values=coles)
    print(s2)
dev.off()

# # bray curtis distance and PCOA
# mx_norm.bray <- vegdist(t(mx_norm), "bray") # bray-curtis distance
# mx_norm.b.pcoa <- cmdscale(mx_norm.bray, k=(nrow(t(mx_norm))-1), eig=TRUE)

# pdf(paste0(res_dir, focus, "_pcoa_abundance_col_by_ds.pdf"), width = 9, height = 7 ) 
# ordiplot(scores(mx_norm.b.pcoa, choices=c(1,2)), type="none", main="PCoA (bray-curtis distance)")
# points(scores(mx_norm.b.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_ds[rownames(scores(mx_norm.b.pcoa, choices=c(1,2)))])
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# legend("topleft", legend=names(colds), pch=16, col=colds)
# dev.off()

# pdf(paste0(res_dir, focus, "_pcoa_abundance_col_by_st.pdf"), width = 9, height = 7 ) 
# ordiplot(scores(mx_norm.b.pcoa, choices=c(1,2)), type="none", main="PCoA (bray-curtis distance)")
# points(scores(mx_norm.b.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_st[rownames(scores(mx_norm.b.pcoa, choices=c(1,2)))])
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# legend("topleft", legend=names(coles[unique(status_l1)]), pch=16, col=coles[unique(status_l1)])
# dev.off()


# compo.bray <- vegdist(compo, "bray")
# compo.bray.pcoa <- cmdscale(compo.bray, k=(nrow(compo)-1), eig=TRUE)
# pdf(paste0(res_dir, focus, "_pcoa_compo_col_by_ds.pdf"), width = 7, height = 7 ) 
# ordiplot(scores(compo.bray.pcoa, choices=c(1,2)), type="none", main="PCoA (bray-curtis distance)")
# points(scores(compo.bray.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_ds[rownames(scores(compo.bray.pcoa, choices=c(1,2)))])
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# legend("topleft", legend=names(colds), pch=16, col=colds)
# dev.off()

# pdf(paste0(res_dir, focus, "_pcoa_compo_col_by_st.pdf"), width = 9, height = 7 ) 
# ordiplot(scores(compo.bray.pcoa, choices=c(1,2)), type="none", main="PCoA (bray-curtis distance)")
# points(scores(compo.bray.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_st[rownames(scores(compo.bray.pcoa, choices=c(1,2)))])
# abline(h=0, lty=3)
# abline(v=0, lty=3)
# legend("topleft", legend=names(coles[unique(status_l1)]), pch=16, col=coles[unique(status_l1)])
# dev.off()


# KO by genomes
kbg <- as.matrix(read.table("ko_by_genome.tsv", sep="\t", header = TRUE))
# not all mgs found in our samples are in the annotation matrix
kbg = kbg[colnames(kbg) %in% rownames(mx_norm),]
# not all mgs found in the cazy-subtrate matrix are in our samples
abundance = mx_norm[rownames(mx_norm) %in% colnames(kbg),]
kbg = kbg[,rownames(abundance)]
# test: there should be only TRUE values when testing the following equality (order of rows)
table(colnames(kbg) == rownames(abundance))
# multiply the matrices to get the KOs for each sample
kbs = kbg %*% abundance 
# bray curtis distance then PCoA
kbs.bray <- vegdist(t(kbs), "bray")
kbs.bray.pcoa <- cmdscale(kbs.bray, k=(nrow(t(kbs))-1), eig=TRUE)
pdf(paste0(res_dir, focus, "_pcoa_ko_col_by_ds.pdf"), width = 7, height = 7 ) 
ordiplot(scores(kbs.bray.pcoa, choices=c(1,2)), type="none", main="")
points(scores(kbs.bray.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_ds[rownames(scores(kbs.bray.pcoa, choices=c(1,2)))])
abline(h=0, lty=3)
abline(v=0, lty=3)
legend("topright", legend=names(colds), pch=16, col=colds)
dev.off()

pdf(paste0(res_dir, focus, "_pcoa_ko_col_by_st.pdf"), width = 7, height = 7 ) 
ordiplot(scores(kbs.bray.pcoa, choices=c(1,2)), type="none", main="")
points(scores(kbs.bray.pcoa, choices=c(1,2)), pch=16, cex=0.6, col = col_by_st[rownames(scores(kbs.bray.pcoa, choices=c(1,2)))])
abline(h=0, lty=3)
abline(v=0, lty=3)
legend("topright", legend=names(coles[unique(status_l1)]), pch=16, col=coles[unique(status_l1)])
dev.off()
