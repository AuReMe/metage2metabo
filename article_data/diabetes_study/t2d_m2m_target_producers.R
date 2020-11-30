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
study = "mgs" # can be "mgs" or "agora"
############################

metadata_file = 'metadata.csv'
if (study == "mgs"){
    target_producers_file = 'butyrate_prod/M_BUTYRIC_ACID_c_producers_dict.tsv'
    target_producers_phylum_file = 'butyrate_prod/M_BUTYRIC_ACID_c_producers_dict_phylum.tsv'
    tax_file = 'mgs_specI.tax'
    tax = as.data.frame(read.table(tax_file, header=FALSE, sep="\t"))
    names(tax) = c("genome","K","P","C","O","F","G","S")
} else if (study == "agora"){
    target_producers_file = 'butyrate_prod/M_but__91__c__93___producers_dict.tsv'
    target_producers_phylum_file = 'butyrate_prod/M_but__91__c__93___producers_dict_phylum.tsv'
    tax_file = 'butyrate_prod/taxonomy_agora.tsv'
    tax = as.data.frame(read.table(tax_file, header=FALSE, sep="\t"))
    names(tax) = c("genome","K","P","C","O","F","G","S")
}

res_dir = paste0('figures_results_', focus, "_", study, "_target_producers_but", "/")
dir.create(file.path(res_dir), showWarnings = FALSE)

# we want to drop samples that had less than 1000 reads passing quality checks
samples_to_be_dropped = c("NG5425286", "NG5425485", "MH0035", "MH0037", "MH0038", "MH0047", "MH0054", "MH0057", "MH0058", "MH0064", "MH0065", "MH0067", "MH0073", "MH0077", "MH0082", "MH0084")

# read input files
target_producers = as.data.frame(read.table(target_producers_file, header=TRUE, sep="\t"))
target_producers_phylum = as.data.frame(read.table(target_producers_phylum_file, header=TRUE, sep="\t"))

# read metadata
metadata = as.data.frame(read.table(metadata_file, header=TRUE, sep=","))

# remove samples that we dropped in each file
metadata = metadata[! metadata$sample %in% samples_to_be_dropped, ]
target_producers = target_producers[! target_producers$sample %in% samples_to_be_dropped, ]
target_producers_phylum = target_producers_phylum[! target_producers_phylum$sample %in% samples_to_be_dropped, ]


mhd_samples = metadata[metadata$subset == 'MHD', 'sample']
swe_samples = metadata[metadata$subset == 'SWE', 'sample']

if (focus=="mhd"){
    # drop the SWE samples
    metadata = metadata[! metadata$sample %in% swe_samples,]
    target_producers = target_producers[! target_producers$sample %in% swe_samples,]
    target_producers_phylum = target_producers_phylum[! target_producers_phylum$sample %in% swe_samples,]
} else if (focus=="swe"){
    # drop the MHD samples
    metadata = metadata[! metadata$sample %in% mhd_samples,]
    target_producers = target_producers[! target_producers$sample %in% mhd_samples,]
    target_producers_phylum = target_producers_phylum[! target_producers_phylum$sample %in% mhd_samples,]
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


########################## Data treatment
# update factors levels in samples
target_producers_phylum$sample <- droplevels(target_producers_phylum$sample)
target_producers$sample <- droplevels(target_producers$sample)

# merge all Firmicutes
target_producers_phylum[startsWith(as.character(target_producers_phylum$producer), "Firmicutes"), "producer"] <- "Firmicutes"
target_producers_phylum$sample <- droplevels(target_producers_phylum$sample)

target_producers_mx = acast(target_producers, producer~sample, value.var="occurrence", fill=0)
target_producers_phylum_mx = acast(target_producers_phylum, producer~sample, value.var="occurrence", fill=0, fun.aggregate = sum)

tax[tax$genome %in% target_producers$producer & startsWith(as.character(tax$P), "Firmicutes"),]

# analyses of all butyrate producers, taking into account the occurrence for each
but_prod = as.data.frame(rowSums(target_producers_mx))
names(but_prod) = c("number")
but_prod$genome = rownames(but_prod)
but_prod$phylum = tax[match(but_prod$genome, tax$genome), "P"]
but_prod[startsWith(as.character(but_prod$phylum), "Firmicutes"), "phylum"] <- "Firmicutes"
group_by(but_prod, phylum) %>%
  summarise(
    count = n(),
    mean = mean(number, na.rm = TRUE),
    sd = sd(number, na.rm = TRUE),
    percent = sum(number)/sum(but_prod$number)*100
  )
rowSums(target_producers_phylum_mx)/ sum(rowSums(target_producers_phylum_mx))*100



# number of members of each phylum among unique tarhet producers
summary(tax[tax$genome %in% unique(target_producers$producer), "P"])

# global diversity of all taxa appearing in at least one community
summary(tax$P)/sum(summary(tax$P))*100
full_diversity = summary(tax$P)
# full_diversity
full_diversity_pc = summary(tax$P)/sum(summary(tax$P))*100
pie(full_diversity_pc, labels = names(full_diversity))
full_diversity = as.data.frame(full_diversity)
names(full_diversity) = c("occurrence")

# full_diversity <- full_diversity %>%
#   arrange(desc(phylum)) %>%
#   mutate(lab.ypos = cumsum(occurrence) - 0.5*occurrence)

# combine Firmicutes rows for mgs
if (study == "mgs"){
    full_diversity = t(full_diversity)
    full_diversity[,"Firmicutes"] = full_diversity[,"Firmicutes_A"] + full_diversity[,"Firmicutes_C"] + full_diversity[,"Firmicutes"]
    full_diversity = full_diversity[, ! colnames(full_diversity) %in% c("Firmicutes_C","Firmicutes_A")]
    full_diversity = as.data.frame(full_diversity)
    names(full_diversity) = c("occurrence")

    # unique targets producers: taxonomy
    utp = summary(tax[tax$genome %in% unique(target_producers$producer), "P"])
    # percentaage of each phylum among all targets producers
    utp = utp/sum(utp)
    utp
    sum(utp[c("Firmicutes", "Firmicutes_A", "Firmicutes_C")])
}

full_diversity$phylum = rownames((full_diversity))

################### Figure
taxcolors = c("#1f77b4ff", "#ff7f0eff", "#ff7f0eff", "#2ca02cff", "#d62728ff", "#d62728ff", "#9467bdff", "#8c564bff", "#e377c2ff", "#e377c2ff", "#7f7f7fff", "#bcbd22ff", "#17becfff", "#57a9e2ff", "#57a9e2ff", "#ffb380ff", "#ffeeaaff", "#bcd35fff", "#ff80e5ff")
names(taxcolors) = c("Firmicutes", "Bacteroidota", "Bacteroidetes", "Proteobacteria", "Actinobacteriota", "Actinobacteria", "Cyanobacteria", "Desulfobacterota_A", "Verrucomicrobiota", "Verrucomicrobia", "Euryarchaeota", "Synergistota", "Spirochaetota", "Spirochaetes", "Thermoplasmatota", "Fusobacteria", "Tenericutes", "Melainabacteria", "Thaumarchaetota")

full_diversity$colors = taxcolors[match(full_diversity$phylum, names(taxcolors))]
full_diversity = full_diversity[order(full_diversity$occurrence, decreasing = TRUE),]
# full_diversity$lab.ypos[full_diversity$occurrence < 14] <- NA

library(plotly)

p <- plot_ly(full_diversity, labels = ~phylum, values = ~occurrence, type = 'pie',textposition = 'outside',textinfo = 'label+value', sort = TRUE, marker = list(colors=full_diversity$colors)) %>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         margin = list(b = 180, l = 50, r = 50, t = 10)) %>% 
         add_trace(showlegend = FALSE)
setwd(res_dir)
namefig = "all_mn_tax.pdf"
orca(p, namefig)


but_producers_by_sample = as.data.frame(colSums(target_producers_mx))
but_producers_by_sample$sample = rownames(but_producers_by_sample)
names(but_producers_by_sample) = c("number", "sample")
but_producers_by_sample$status_l1 = status_l1[but_producers_by_sample$sample]
but_producers_by_sample$status_l2 = status_l2[but_producers_by_sample$sample]
but_producers_by_sample$status_l3 = status_l3[but_producers_by_sample$sample]

# Anova tests effect of status on number of putative butyrate producers
fit1 = aov(number ~ status_l1, but_producers_by_sample)
fit2 = aov(number ~ status_l2, but_producers_by_sample)
fit3 = aov(number ~ status_l3, but_producers_by_sample)

# summary(fit1)
# TukeyHSD(fit1)
summary(fit2)
TukeyHSD(fit2)
etaSquared(fit2)
# summary(fit3)
# TukeyHSD(fit3)

group_by(but_producers_by_sample, status_l1) %>%
  summarise(
    count = n(),
    mean = mean(number, na.rm = TRUE),
    sd = sd(number, na.rm = TRUE)
  )
group_by(but_producers_by_sample, status_l2) %>%
  summarise(
    count = n(),
    mean = mean(number, na.rm = TRUE),
    sd = sd(number, na.rm = TRUE)
  )
group_by(but_producers_by_sample, status_l3) %>%
  summarise(
    count = n(),
    mean = mean(number, na.rm = TRUE),
    sd = sd(number, na.rm = TRUE)
  )
