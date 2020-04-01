###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(compositions)
library(ape)
library(factoextra)
library(RColorBrewer)
library(philr)
library(phyloseq)
library(UpSetR)
library(reshape2)
library(vegan)
library(phylofactor)
library(ggtree)

##########
#ENV SETUP
##########
PATH="/home/lymelab/mann/edgar"
setwd(PATH)
cols <- brewer.pal(6, "Set2")

###########
#LOAD DATA
###########
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "edgar_map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
rownames(rawmetadata) <- rawmetadata$SampleID
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("taxonomy_L7.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID) #record samples absent in either metadata or OTU table
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
tree <- read.tree("rep_set.root.tre") #load representative tree

#prune out Heather's samples
seqtab.filtered <- seqtab.filtered[!row.names(seqtab.filtered) %in% notinmeta,]

#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

################
#PHILR DISTANCE 
################
#philr transform for full dataset
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

# Heirarchical cluster dendrogram
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)

hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$Season))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/philr_dendrogram_season.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# PCA
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))

pdf("figs/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Season, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() 
dev.off()

pdf("figs/philr_pca_lysing_matrix.pdf")
fviz_pca_ind(pca, habillage=merge$Lysing_matrix, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() 
dev.off()

############
#Upset plot
############
map <- as.matrix(read.table("edgar_map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merge$Season), FUN=sum)
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- data.frame(t(agg[,-1]))
#upsetR 
pdf("figs/upset.pdf", onefile=F)
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs per season", mb.ratio = c(0.55, 0.45))
dev.off()


#########
#Adonis 
#########
# is there a difference in microbial diversity across season?
metadata <- as(sample_data(ps.dat), "data.frame")

adonis(philr.dist ~ Season, data=metadata)

# Call:
# adonis(formula = philr.dist ~ Season, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Season     2    3802.9 1901.47  6.6641 0.29404  0.001 ***
# Residuals 32    9130.5  285.33         0.70596
# Total     34   12933.4                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

