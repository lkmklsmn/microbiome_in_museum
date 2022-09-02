# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)

setwd("/Users/lukas.simon/Downloads/lounsbery")

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

# Load data ####
meta <- read_excel("Sampling list_08102020.xlsx", sheet = 1)
meta$id <- paste("L", meta$`Sampling number`, sep = "_")

table <- read_qza("table.qza")
table <- as.matrix(table$data)

# Restrict to samples present in both ####
ok <- intersect(meta$id, colnames(table))
samples <- sample_data(meta[match(ok, meta$id),])
OTU <- otu_table(table[, ok], taxa_are_rows = T)
sample_names(OTU) <- sample_names(samples)

# Load taxonomy info ####
classified <- read_qza("classified_rep_seqs.qza")
tmp <- classified$data
tmp <- as.character(tmp$Taxon)
tmp <- do.call(rbind, lapply(tmp, function(x){
  empty <- rep(NA, 7)
  if(x == "Unassigned") return(empty)
  taxons <- strsplit(x, ";", fixed = T)[[1]]
  empty[1:length(taxons)] <- taxons
  empty
}))
rownames(tmp) <- rownames(OTU)
colnames(tmp) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tmp <- as.matrix(tmp)
TAX <- tax_table(tmp)

# Create phyloseq object ####
ok <- which(!is.na(tmp[, "Phylum"]))

# Create phyloseq object ####
p <- phyloseq_miko(OTU, TAX, samples)

# Remove taxa with NA for phylum ####
ok <- names(which(!is.na(TAX@.Data[, "Phylum"])))
p <- prune_taxa(taxa = ok, p)

# Create some plots ####
plot_bar(p, fill = "Phylum")

tmp <- log(OTU@.Data + 1)
tmp <- tmp[which(apply(tmp, 1, var) > 0),]
pca <- prcomp(t(tmp), scale. = T)


counts <- OTU@.Data

aframe <- data.frame(pca$x, data.frame(samples), total_counts = colSums(counts))
ggplot(aframe, aes(Collection, total_counts, color = Collection)) +
  geom_boxplot() + geom_point() +
  theme_bw()

ggplot(aframe, aes(Object, log10(total_counts), color = Object)) +
  geom_boxplot() + geom_point() +
  theme_bw()

info <- data.frame(samples)
info$total <- colSums(counts)

ok <- which(info$Object %in% c("Ischta Gate, right side", "Procession Street, left side") & info$total > 10000)
ok <- which(info$total > 10000)
ok <- rownames(info)[ok]
calc_sums <- function(counts, taxa = "Phylum"){
  asplit <- split(rownames(counts), TAX@.Data[, taxa])
  sums <- do.call(rbind, lapply(asplit, function(x) colSums(counts[x, ])))
  sums
}
sums <- calc_sums(counts[, ok], taxa = "Domain")

aframe <- data.frame(t(sums), info[ok,])
ggplot(aframe, aes(log10(D_1__Cyanobacteria + 1), log10(D_1__Nanoarchaeaeota + 1), color = spot)) +#, D_1__Proteobacteria))
  geom_point() +
  theme_bw()

Down_Sample_Matrix <- function (expr_mat)
{
min_lib_size <- min(colSums(expr_mat))
down_sample <- function(x) {
prob <- min_lib_size/sum(x)
return(unlist(lapply(x, function(y) {
rbinom(1, y, prob)
})))
}
down_sampled_mat <- apply(expr_mat, 2, down_sample)
return(down_sampled_mat)
}

library(edgeR)

tmp <- log(Down_Sample_Matrix(counts[, ok]) + 1)
tmp <- tmp[which(apply(tmp, 1, var) > 0), ]
pca <- prcomp(t(tmp), scale. = T)
aframe <- data.frame(info[ok,], t(sums))

ok <- which(info$total > 1000)
ok <- rownames(info)[ok]
sums <- calc_sums(counts[, ok], taxa = "Domain")
aframe <- data.frame(info[ok,], t(sums))
ggplot(aframe, aes(D_0__Archaea, D_0__Bacteria, color = Institution)) +
  geom_point() +
  theme_bw()

aframe$touched <- "touched"
aframe$touched[grep("untouched", aframe$comments)] <- "untouched"
ggplot(aframe, aes(D_0__Archaea, D_0__Bacteria, color = Object, shape = touched)) +
  facet_wrap(~ Institution) +
  geom_point() +
  theme_bw()

sums <- calc_sums(counts[, ok], taxa = "Species")
detected <- apply(counts, 1, function(x) sum(x > 0))

prop <- t(t(sums)/colSums(sums))
pca <- prcomp(t(prop))
aframe <- data.frame(info[ok,], pca$x)
ggplot(aframe[aframe$Object %in% c("Procession Street, left side", "Ischta Gate, right side"),],
       aes(PC1, PC2, color = comments, label = comments)) +
  geom_point() + geom_label() +
  theme_bw()

counts[which(!is.na(TAX@.Data[,1])),]

