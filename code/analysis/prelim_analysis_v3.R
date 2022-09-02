# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(microbiome)
library(colorspace)

setwd("/Users/lukas/OneDrive/Documents/GitHub/lounsbery/")

# Define some functions ####
phyloseq_miko <- function (...) 
{
  arglist <- list(...)
  names(arglist) <- NULL
  arglist <- arglist[sapply(arglist, phyloseq:::is.component.class)]
  splatlist <- sapply(arglist, phyloseq:::splat.phyloseq.objects)
  splatlist = lapply(splatlist, function(x) {
    taxa_names(x) <- gsub("\"", "", taxa_names(x), 
                          fixed = TRUE)
    taxa_names(x) <- gsub("'", "", taxa_names(x), 
                          fixed = TRUE)
    return(x)
  })
  if (length(splatlist) > length(phyloseq:::get.component.classes())) {
    stop("Too many components provided\n")
  }
  else if (length(names(splatlist)) > length(unique(names(splatlist)))) {
    stop("Only one of each component type allowed.\n", 
         "For merging multiple objects of the same type/class, try merge_phyloseq(...)\n")
  }
  else if (length(splatlist) == 1) {
    return(arglist[[1]])
  }
  else {
    ps <- do.call("new", c(list(Class = "phyloseq"), 
                           splatlist))
  }
  shared_taxa = phyloseq:::intersect_taxa(ps)
  shared_samples = phyloseq:::intersect_samples(ps)
  if (length(shared_taxa) < 1) {
    stop("Problem with OTU/taxa indices among those you provided.\n", 
         "Check using intersect() and taxa_names()\n")
  }
  if (length(shared_samples) < 1) {
    stop("Problem with sample indices among those you provided.\n", 
         "Check using intersect() and taxa_names()\n")
  }
  ps = prune_taxa(shared_taxa, ps)
  #ps = prune_samples(shared_samples, ps)
  ps = phyloseq:::index_reorder(ps, "both")
  if (!is.null(phy_tree(ps, FALSE))) {
    ps@phy_tree <- fix_phylo(ps@phy_tree)
  }
  return(ps)
}

# Load sample sheet ####
meta <- read.delim("sample_sheet.tsv")
meta$object <- gsub("Mobbi Collection, Müggelsee 1940", "Mobbi collection", meta$object)
meta <- meta[!is.na(meta$sample.id),]
info <- sample_data(meta)

# Load OTU counts ####
table <- read_qza("data/table.qza")
table <- as.matrix(table$data)
OTU <- otu_table(table, taxa_are_rows = T)
sample_names(info) <- sample_names(OTU)

# Load taxonomy info ####
classified <- read_qza("data/taxonomy_downloaded.qza")
tmp <- classified$data
tmp <- as.character(tmp$Taxon)
tmp <- do.call(rbind, lapply(tmp, function(x){
  empty <- rep(NA, 7)
  if(x == "Unassigned") return(empty)
  taxons <- strsplit(x, "; ", fixed = T)[[1]]
  empty[1:length(taxons)] <- taxons
  empty
}))
rownames(tmp) <- classified$data$Feature.ID
colnames(tmp) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tmp <- as.matrix(tmp)
TAX <- tax_table(tmp)

# Do dim reduction manually ####
counts <- as.matrix(read_qza("data/table.qza")$data)
total <- colSums(counts)
num_detected <- apply(counts, 2, function(x) sum(x > 0))

aframe <- data.frame(total, num_detected, meta)
ggplot(aframe, aes(total + 1, num_detected,
                   color = weight, shape = institution,
                   label = sample.id)) +
  geom_point() + geom_label() +
  scale_x_continuous(trans = "log10") +
  ylab("Number of feature detected") +
  xlab("Total number of reads") +
  theme_bw()

zeros <- apply(counts, 1, function(x) sum(x > 0))
counts <- counts[which(zeros > 1),]

normalized <- t(t(counts)/colSums(counts))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe <- data.frame(meta, total = colSums(counts), pca$x[, 1:10])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

p1 <- ggplot(aframe, aes(PC1, PC2, color = institution, shape = institution, label = sample.id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

p2 <- ggplot(aframe, aes(PC1, PC2, color = weight, shape = institution, label = sample.id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

grid.arrange(p1, p2, ncol = 2)


ok <- which(meta$weight == "ok")
pca <- prcomp(t(normalized[, ok]))
aframe <- data.frame(meta[ok, ], total = colSums(counts)[ok], pca$x[, 1:10])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

ggplot(aframe, aes(PC1, PC2, color = institution, shape = institution, label = sample.id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

# HMP ####
temp_otutax = tempfile()
otu_url = "http://downloads.hmpdacc.org/data/HMQCP/otu_table_psn_v35.txt.gz"
download.file(otu_url, temp_otutax)

# Create phylum barplot ####
level <- "Phylum"
ok <- which(colSums(OTU) > 1000)
tmp <- t(t(OTU[, ok])/colSums(OTU[, ok]))
asplit <- split(rownames(TAX), TAX[, level])
norm <- do.call(rbind, lapply(asplit, function(x) colSums(OTU[x, ok])))
norm <- t(t(norm)/colSums(norm))

if(level=="Genus") top_hits <- tail(setdiff(names(sort(rowMeans(norm))), "g__"), 15)
if(level=="Phylum") top_hits <- tail(setdiff(names(sort(rowMeans(norm))), "p__"), 10)
rest <- colSums(norm[setdiff(rownames(norm), top_hits),])

aframe <- data.frame(t(norm[top_hits, ]), rest, meta[ok, ])
ord <- order(aframe$institution, aframe$object)
aframe$sample.id <- factor(aframe$sample.id, levels = (aframe$sample.id[ord]))

aframe <- melt(aframe, id.vars = colnames(meta))

farben <- qualitative_hcl(length(top_hits) + 1, "Dark2")
farben[length(farben)] <- "grey"

p_phylum <- ggplot(aframe[aframe$institution != "DSMZ",],
                   aes(x = sample.id, y = value, fill = variable)) + 
  facet_wrap(~ institution, scales = "free_x") +
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_manual(values = farben, name = "Phylum") +
  ylab("Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
grid.arrange(p1, p_phylum, ncol = 2, widths = c(2, 3))

# Combine into one analysis ####
tmp <- meta[meta$institution == "MfN" & meta$weight == "ok",]
tmp <- meta[meta$institution == "MfN",]
tmp <- tmp[tmp$comments != "reference",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), pca$x)

aframe_tmp$touched <- "yes"
aframe_tmp$touched[grep("untouched", aframe_tmp$comments)] <- "no"
aframe_tmp$groups <- "molluscs"
aframe_tmp$groups[aframe_tmp$collection == "VertPal"] <- "VertPal"
aframe_tmp$groups[aframe_tmp$touched == "yes"] <- "touched"

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

p_pca <- ggplot(aframe_tmp, aes(PC1, PC2, group = groups, color = groups, shape = collection)) +
  geom_point() +
  stat_ellipse(linetype = 'dashed') +
  ggtitle("PCA - MfN") +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  theme_classic()

res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x, aframe_tmp[, c("touched", "object", "collection")])
  afit <- try(summary(aov(expr ~ collection + touched, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  as.numeric(afit[[1]][2, 4:5])
}))
colnames(res) <- c("F_value", "p_value")
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res[, 2]),]

anno_col <- data.frame(aframe_tmp[, c("touched", "collection", "object", "comments")])
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(TAX@.Data[rownames(normalized),])
rownames(anno_row) <- rownames(normalized)

sig <- rownames(res)[p.adjust(res[,2], method = "BH") < 0.1]

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollucs` = "blue", `VertPal` = "darkblue")
object <- c(`Chiton` = "bisque3", `Tendaguru` = "darksalmon", `Triton horn`= "yellow", `Mobbi collection`= "brown")
ann_colors <- list(touched = touched)#,
#                   collection = collection,
#                   object = object)

pheatmap(normalized[sig, order(anno_col$collection, anno_col$touched, anno_col$object)],
         cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = 10,
         annotation_col = anno_col[, c("object", "collection", "touched")],
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")],
         annotation_colors = ann_colors)

# Assess loss of identity ####
res_molluscs_vs_fossils <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x, aframe_tmp[, c("touched", "object", "collection")])
  aframe <- aframe[aframe$touched == "no", ]
  afit <- try(summary(lm(expr ~ collection, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  coefficients(afit)[2, c(1, 4)]
}))
colnames(res_molluscs_vs_fossils) <- c("coef", "p_value")
res_molluscs_vs_fossils <- data.frame(res_molluscs_vs_fossils)
res_molluscs_vs_fossils$padj <- p.adjust(res_molluscs_vs_fossils$p_value, method = "BH")
res_molluscs_vs_fossils <- data.frame(res_molluscs_vs_fossils, TAX@.Data[rownames(res_molluscs_vs_fossils),])
res_molluscs_vs_fossils <- res_molluscs_vs_fossils[sort.list(res_molluscs_vs_fossils[, 2]),]

ggplot(res_molluscs_vs_fossils, aes(coef, -log10(p_value))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Abundance molluscs/fossils") +
  theme_bw()

rownames(aframe_tmp) <- colnames(normalized)

pheatmap(normalized[rownames(res_molluscs_vs_fossils)[1:20], order(aframe_tmp$groups, aframe_tmp$collection)],
         cluster_cols = F,
         gaps_col = c(6, 13),
         annotation_col = aframe_tmp[, c("groups", "collection")],
         show_rownames = F)

aframe_tmp$pc_1 <- NA
aframe_tmp$pc_2 <- NA

ok <- which(aframe_tmp$touched == "no")
tmp <- normalized[, ok]
pca <- prcomp(t(tmp))
tmp <- data.frame(pca$x, aframe_tmp[ok,])
ggplot(tmp, aes(PC1, PC2, color = collection)) +
  geom_point() +
  theme_bw()

aframe_tmp$pc_1[ok] <- pca$x[,1]
aframe_tmp$pc_2[ok] <- pca$x[,2]

ok <- aframe_tmp$touched == "yes"
preds <- predict(pca, t(normalized[, ok]))

aframe_tmp$pc_1[ok] <- preds[,1]
aframe_tmp$pc_2[ok] <- preds[,2]

p_pca_2 <- ggplot(aframe_tmp, aes(pc_1, pc_2, shape = collection, color = groups)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  ggtitle("Principal component analysis") +
  theme_bw()

p_pc1 <- ggplot(aframe_tmp[aframe_tmp$touched == "yes",],
       aes(collection, pc_1, color = collection)) +
  geom_boxplot() + geom_point() +
  ylab("Principal component 1") +
  theme_bw()

grid.arrange(p_pca_2, p_pc1, ncol = 2)


aframe_tmp$expr <- normalized["1d5880b640928ac56a26ab9e58c0ce9a", ]
ggplot(aframe_tmp, aes(groups, expr, color = groups)) +
  geom_boxplot() + geom_point() +
  theme_bw()

plot_feature <- function(feature, title){
  aframe <- data.frame(expr = normalized[feature, ], aframe_tmp[, c("touched", "object", "collection")])
  pval <- summary(aov(expr ~ collection  + touched, data=aframe))[[1]][2,5]
  nom <- TAX@.Data[feature, -1]
  nom <- paste(nom, collapse = ",")
  ggplot(aframe, aes(touched, expr, color = touched)) +
    facet_wrap(~ collection) +
    geom_boxplot() + geom_point() +
    ggtitle(paste(title, "P", signif(pval, 2))) +
    ylab("Abundance") +
    theme_bw()
}
p1 <- plot_feature("0664a0c96e567c98809225119f8fad99", "Proprionibacterium")

# Validation of touch signature in Gate lions ####
tmp <- data.frame(feature = colSums(counts[sig,]),
                  meta,
                  total)

tmp$touched <- "no"
tmp$touched[tmp$comments == "touched"] <- "yes"
tmp <- tmp[grep("Lion", tmp$object),]
tmp <- tmp[tmp$weight == "ok",]

summary(lm(log10(feature + 1) ~ object + touched, data = tmp))

ggplot(tmp, aes(touched, log10(feature + 1), color = touched)) +
  facet_wrap(~ object) +
  geom_boxplot() + geom_point() +
  ggtitle("Gate Lion") +
  ylab("'Touch' signature") +
  theme_bw()

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Streptococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Streptococcus abundance")

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Staphylococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Staphylococcus abundance")


# Combined analysis ####
tmp <- counts[sig, rev(names(sort(colSums(counts[sig,])/colSums(counts))))]
anno_col <- meta
rownames(anno_col) <- meta$sample.id
anno_col$touched <- "no"
anno_col$touched[which(anno_col$comments %in% c("touched", "touched by users", "very much touched by users"))] <- "yes"
pheatmap(log10(tmp + 1), cluster_cols = F,
         annotation_row = anno_row[, c("Genus", "Phylum")], show_rownames = F,
         annotation_col = anno_col[, c('institution', 'collection', 'object', 'touched')])

