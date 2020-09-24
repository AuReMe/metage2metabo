library(ggplot2)
library(ggfortify)
library(plyr)
library(viridis)
iscope_file_gut = "gut_iscope_sizes.tsv" # individual scopes of GSMNs reconstructed from the genomes of cultured gut species
iscope_file_rumen = "rumen_iscope_sizes.tsv" # individual scopes of GSMNs reconstructed from the rumen MAGs
cpd_prod_file_gut = "gut_cpd_producers.tsv" # number of individual producers for each metabolites of the rumen GSMNs
cpd_prod_file_rumen = "rumen_cpd_producers.tsv" # number of individual producers for each metabolites of the (cultured) gut GSMNs
nb_genes_gut_file = "gut_genes_per_genomes.tsv" # number of genes in genomes of cultured gut species
nb_genes_rumen_file = "degradation_impact_genes.tsv" # number of genes in rumen MAGs with or without degradation
t2d_mgs_file = "reconstats_mgs.tsv" # statistics on GSMNs reconstructed in the Diabetes dataset (MGS)
t2d_speci_file = "reconstats_specI.tsv" # statistics on GSMNs reconstructed in the Diabetes dataset (from SpecI genomes)

iscope_gut = as.data.frame(read.table(iscope_file_gut, header=TRUE, sep="\t", check.names=FALSE))
iscope_rumen = as.data.frame(read.table(iscope_file_rumen, header=TRUE, sep="\t", check.names=FALSE))
cpd_prod_gut = as.data.frame(read.table(cpd_prod_file_gut, header=TRUE, sep="\t", check.names=FALSE))
cpd_prod_rumen = as.data.frame(read.table(cpd_prod_file_rumen, header=TRUE, sep="\t", check.names=FALSE))
allgenes_gut = as.data.frame(read.table(nb_genes_gut_file, header=TRUE, sep="\t", check.names=FALSE))
allgenes_rumen = as.data.frame(read.table(nb_genes_rumen_file, header=TRUE, sep="\t", check.names=FALSE))
t2d_specI = as.data.frame(read.table(t2d_speci_file, header=TRUE, sep="\t", check.names=FALSE))
t2d_mgs = as.data.frame(read.table(t2d_mgs_file, header=TRUE, sep="\t", check.names=FALSE))

t2s_stats = rbind(t2d_mgs, t2d_specI)

allgenes_rumen = allgenes_rumen[,c("genome", "original")]
names(allgenes_rumen) = c("genome", "number_genes")

res_dir = paste0('supp_fig_r', "/")
dir.create(file.path(res_dir), showWarnings = FALSE)


# cpd_prod_gut = cpd_prod_gut[order(cpd_prod_gut$producers_nb, decreasing=TRUE),]
cpd_prod_gut_plot = ggplot(data=cpd_prod_gut, aes(x=reorder(compounds, -producers_nb), y=producers_nb)) +geom_bar(stat = "identity",color=viridis(3)[1], fill=viridis(3)[1]) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.x = element_blank(),
axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
legend.text=element_text(size=24)) + scale_y_continuous(expand = c(0, 0)) +
labs(x = paste0("metabolites (", length(cpd_prod_gut$compounds), ")"), y = "Number of GSMN producing the metabolite")

pdf(paste0(res_dir, "cpd_scope_nb_gut.pdf"), width = 11, height = 7)
cpd_prod_gut_plot
dev.off()

# cpd_prod_gut = cpd_prod_gut[order(cpd_prod_gut$producers_nb, decreasing=TRUE),]
cpd_prod_rumen_plot = ggplot(data=cpd_prod_rumen, aes(x=reorder(compounds, -producers_nb), y=producers_nb)) +geom_bar(stat = "identity",color=viridis(3)[2], fill=viridis(3)[2]) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.x = element_blank(),
axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=26), legend.position="none",
legend.text=element_text(size=24)) + scale_y_continuous(expand = c(0, 0)) +
labs(x = paste0("metabolites (", length(cpd_prod_rumen$compounds), ")"), y = "Number of GSMN producing the metabolite")

pdf(paste0(res_dir, "cpd_scope_nb_rumen.pdf"), width = 11, height = 7)
cpd_prod_rumen_plot
dev.off()

iscope_gut_plot = ggplot(data=iscope_gut, aes(x=reorder(bacterium, -scope_size), y=scope_size)) +geom_bar(stat = "identity",color=viridis(3)[1], fill=viridis(3)[1]) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.x = element_blank(),
axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
legend.text=element_text(size=24)) + scale_y_continuous(expand = c(0, 0)) +
labs(x = paste0("bacteria (", length(iscope_gut$bacterium), ")"), y = "Size of individual metabolic potential")+
geom_hline(aes(yintercept=93), color="white", linetype="dashed")

pdf(paste0(res_dir, "size_scope_nb_gut.pdf"), width = 11, height = 7)
iscope_gut_plot
dev.off()

iscope_rumen_plot = ggplot(data=iscope_rumen, aes(x=reorder(bacterium, -scope_size), y=scope_size)) +geom_bar(stat = "identity",color=viridis(3)[2], fill=viridis(3)[2]) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.x = element_blank(), legend.position="none",
axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), 
legend.text=element_text(size=24)) + scale_y_continuous(expand = c(0, 0)) +
labs(x = paste0("bacteria (", length(iscope_rumen$bacterium), ")"), y = "Size of individual metabolic potential")+
geom_hline(aes(yintercept=26), color="white", linetype="dashed")

pdf(paste0(res_dir, "size_scope_nb_rumen.pdf"), width = 11, height = 7)
iscope_rumen_plot
dev.off()

##################################################################################################################

padmetstats_gut_file <- "padmet_stats_gut.tsv" # statistics on GSMNs reconstructed from the genomes of cultured gut species
pstats_gut = as.data.frame(read.table(padmetstats_gut_file, header=TRUE, sep="\t", row.names = 1))
pstats_gut$origin = "gut"

padmetstats_rumen_file <- "padmet_stats_rumen.tsv" # statistics on GSMNs reconstructed from the rumen MAGs
pstats_rumen = as.data.frame(read.table(padmetstats_rumen_file, header=TRUE, sep="\t", row.names = 1))
pstats_rumen$origin = "rumen"

# concatenate both dataframes
pstats_all = rbind(pstats_gut, pstats_rumen)

# means by condition
cdat <- ddply(pstats_all, "origin", summarise, pwy.mean=mean(pathways), rxn.mean=mean(reactions), rxngene.mean=mean(reactions_with_gene_association), cpd.mean=mean(compounds), genes.mean=mean(genes))

allgenes_rumen$origin = "rumen"
allgenes_gut$origin = "gut"
nb_genes_all = rbind(allgenes_rumen, allgenes_gut)

p = ggplot(pstats_all, aes(x=pathways, fill=origin, color=origin)) +
    geom_histogram(binwidth=5, alpha=.7, position="identity" ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
    legend.text=element_text(size=24)) +
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) + 
    scale_y_sqrt(breaks=c(1,10,30,50,100))
    #geom_vline(data=cdat, aes(xintercept=pwy.mean),
    # linetype="dashed", size=1, colour="black") +

pdf(paste0(res_dir, "pathways.pdf"), width = 11, height = 7)
p
dev.off()

q = ggplot(pstats_all, aes(x=reactions, fill=origin, color=origin)) +
    geom_histogram(binwidth=30, alpha=.7, position="identity") +
    # geom_vline(data=cdat, aes(xintercept=rxn.mean),
    #            linetype="dashed", size=1, colour="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
    legend.text=element_text(size=24))+
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) 

pdf(paste0(res_dir, "reactions.pdf"), width = 11, height = 7)
q
dev.off()

r = ggplot(pstats_all, aes(x=compounds, fill=origin, color=origin)) +
    geom_histogram(binwidth=30, alpha=.7, position="identity") +
    # geom_vline(data=cdat, aes(xintercept=cpd.mean),
    #            linetype="dashed", size=1, colour="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
    legend.text=element_text(size=24))+
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) 

pdf(paste0(res_dir, "compounds.pdf"), width = 11, height = 7)
r
dev.off()

s = ggplot(pstats_all, aes(x=reactions_with_gene_association, fill=origin, color=origin)) +
    geom_histogram(binwidth=30, alpha=.7, position="identity") +
    # geom_vline(data=cdat, aes(xintercept=rxngene.mean),
    #            linetype="dashed", size=1, colour="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.position="none",
    legend.text=element_text(size=24))+
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) +
    labs(x = "reactions with gene-associations")

pdf(paste0(res_dir, "reaction_with_genes.pdf"), width = 11, height = 7)
s
dev.off()

u = ggplot(pstats_all, aes(x=genes, fill=origin, color=origin)) +
    geom_histogram(binwidth=20, alpha=.7, position="identity") +
    # geom_vline(data=cdat, aes(xintercept=rxngene.mean),
    #            linetype="dashed", size=1, colour="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16),  legend.position="none",
    legend.text=element_text(size=24)) +
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) +
    labs(x = "genes in the GSMN")

pdf(paste0(res_dir, "genes_gsmn.pdf"), width = 11, height = 7)
u
dev.off()

# genes in the genomes
v = ggplot(nb_genes_all, aes(x=number_genes, fill=origin, color=origin)) +
    geom_histogram(binwidth=100, alpha=.7, position="identity") +
    # geom_vline(data=cdat, aes(xintercept=rxngene.mean),
    #            linetype="dashed", size=1, colour="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
     axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16),  legend.position="none",
    legend.text=element_text(size=24)) +
    scale_color_manual(values=viridis(3)[1:2]) +
    scale_fill_manual(values=viridis(3)[1:2]) +
    labs(x = "genes in the genome")

pdf(paste0(res_dir, "genes_genome.pdf"), width = 11, height = 7)
v
dev.off()


pstats_norm = pstats_all
library(scales)
pstats_norm = pstats_all
pstats_norm$pathways <- rescale(pstats_norm$pathways)
pstats_norm$compounds <- rescale(pstats_norm$compounds)
pstats_norm$reactions <- rescale(pstats_norm$reactions)
pstats_norm$reactions_with_gene_association <- rescale(pstats_norm$reactions_with_gene_association)

df <- pstats_all[c(1,2,3,4,5)]
pca_plot <- autoplot(prcomp(df), data=pstats_all, colour="origin") + scale_color_manual(values=viridis(3)[1:2]) + theme_bw() + theme(axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), 
    legend.text=element_text(size=24),legend.position="none")#,
        #  loadings = TRUE, loadings.colour = 'blue',
        #  loadings.label = TRUE, loadings.label.size = 6, label.color = "black") 

pdf(paste0(res_dir, "pca.pdf"), width = 11, height = 7)
pca_plot
dev.off()

pca_plot <- autoplot(prcomp(df), data=pstats_all, colour="origin", loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 6, label.color = "black") + scale_color_manual(values=viridis(3)[1:2]) + theme_bw() + theme(axis.text=element_text(size=24), axis.title=element_text(size=26),legend.title=element_text(size=16), legend.text=element_text(size=24)) 
