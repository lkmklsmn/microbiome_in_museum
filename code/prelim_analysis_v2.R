# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)

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

# Create phyloseq object ####
p <- phyloseq_miko(OTU, TAX, info)

# Create UNIFRAC plot ####
ok = genefilter_sample(p, filterfun_sample(function(x) x > 5), A=0.05*nsamples(p))
GP1 = prune_taxa(ok, p)

ok = is.na(TAX[,"Phylum"])
GP1 = prune_taxa(ok, GP1)

ok = names(which(colSums(OTU) > 100))
GP1 = prune_samples(ok, GP1)

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(physeq = GP1, method = "MDS")
plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
plot_ordination(GP1, GP.ord, type = "samples", color = "Object.number", shape = 'Institution', title = "samples") 

# Plot PCA ####
aframe <- data.frame(GP.ord$vectors, sample_data(GP1))
aframe$label <- NA
aframe$label[aframe$type.of.room == ""] <- aframe$sample.id[aframe$type.of.room == ""]

ggplot(aframe, aes(Axis.1, Axis.2, shape = Institution, color = Object, label = label)) +
  geom_point() + ggrepel::geom_label_repel() +
  theme_bw()

ggplot(aframe, aes(Axis.1, Axis.2, shape = Institution, color = spot == "Air")) +
  geom_point() +
  theme_bw()

ggplot(aframe, aes(Axis.3, Axis.4, shape = Institution, color = factor(Object.number))) +
  geom_point() +
  theme_bw()

# Do dim reduction manually ####
counts <- as.matrix(read_qza("data/table.qza")$data)

zeros <- apply(counts, 1, function(x) sum(x > 0))
counts <- counts[which(zeros > 1),]

normalized <- t(t(counts)/colSums(counts))
normalized <- sqrt(normalized)

#library(edgeR)
#normalized <- log(cpm(counts) + 1)

#normalized <- t(t(counts + 1)/colSums(counts))*1e6
#normalized <- log(normalized)

pca <- prcomp(t(normalized))#, scale. = T)
aframe <- data.frame(meta, total = colSums(counts), pca$x[, 1:10])

mds <- cmdscale(dist(t(normalized)), k = 2)
aframe <- data.frame(meta, total = colSums(counts), pca$x[, 1:10], mds)

ggplot(aframe, aes(X1, X2, color = Institution, shape = Institution, label = sample.id)) +
  geom_point() + geom_label() +
  theme_classic()

ggplot(aframe, aes(X1, X2, color = log10(total), shape = Institution)) +
  geom_point() +
  theme_classic()

ggplot(aframe, aes(PC1, PC2, color = log10(total), shape = Institution)) +
  geom_point() +
  theme_classic()

ggplot(aframe, aes(PC1, PC2, color = Institution, shape = Institution, label = sample.id)) +
  geom_point() + geom_label() +
  theme_classic()

# Run positive unlabeled learning ####
outcome <- meta$Institution == "Control"
positives <- which(outcome)
unlabeled <- which(!outcome)
predictions <- do.call(rbind, lapply(1:100, function(x){
  train <- c(positives, sample(unlabeled, length(positives)))
  #afit <- lda(x = t(normalized[, train]), grouping = outcome[train])  
  aframe <- data.frame(outcome = factor(outcome), t(normalized))
  rf <- randomForest::randomForest(outcome ~., data = aframe[train, ])
  preds <- predict(rf, aframe)
  preds <- as.numeric(preds)
  preds[setdiff(train, positives)] <- NA
  preds
}))

avg_pred <- apply(predictions, 2, function(x) mean(x, na.rm = T))
meta$avg_pred <- avg_pred
aframe <- data.frame(avg_pred, total = colSums(counts), meta)

p1 <- ggplot(aframe, aes(avg_pred, log10(total), color = Institution == "Control", shape = Institution)) +
  geom_point() +
  ylab("Total count [log10]") + xlab("Average prediction") +
  ggtitle("Positive-unlabeled learning") +
  theme_bw()

p2 <- ggplot(aframe, aes(avg_pred, log10(total), color = spot == "Air", shape = Institution)) +
  geom_point() +
  ylab("Total count [log10]") + xlab("Average prediction") +
  ggtitle("Positive-unlabeled learning") +
  theme_bw()

gridExtra::grid.arrange(p1, p2, ncol = 2)

# Run analysis on good lion samples ####
ok <- which(aframe$avg_pred < 1.5)
tmp <- meta[ok,]
ok_samples <- tmp$sample.id[tmp$Object %in% c("Lion left", "Lion right")]
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
aframe_tmp <- data.frame(meta[match(ok_samples, meta$sample.id), ], total = colSums(counts_tmp), mds)

ggplot(aframe_tmp, aes(X1, X2, color = touched, shape = Object)) +
  geom_point() +
  ggtitle("MDS - Lion") +
  theme_classic()

# Hellenistic Haal ####
ok_samples <- tmp$sample.id[tmp$Object %in% c("Hellenistic Haal")]
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
aframe_tmp <- data.frame(meta[match(ok_samples, meta$sample.id), ], total = colSums(counts_tmp), mds)

ggplot(aframe_tmp, aes(X1, X2, label = spot)) +
  geom_point() + ggrepel::geom_label_repel() +
  ggtitle("MDS - Hellenistic Haal") +
  theme_classic()

# Tendaguru ####
tmp <- meta[which(meta$Object == "Tendaguru"), ]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(ok_samples, meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = comments)) +
  geom_point() +
  ggtitle("MDS - Tendaguru") +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = comments)) +
  geom_point() +
  ggtitle("PCA - Tendaguru") +
  theme_classic()

aframe_tmp <- data.frame(PC2 = pca$rotation[,2], TAX@.Data[rownames(pca$rotation), ])
asplit <- split(aframe_tmp$PC2, aframe_tmp$Genus)
aframe_tmp$Genus <- factor(aframe_tmp$Genus, levels = names(sort(unlist(lapply(asplit, median)))))
ggplot(aframe_tmp[!is.na(aframe_tmp$Genus),], aes(PC2, Genus, color = Phylum)) +
  geom_boxplot() +
  xlab("PC2 loading") +
  theme_bw()


# Run diff expr
library(DESeq2)

des <- DESeqDataSetFromMatrix(countData = Down_Sample_Matrix(counts_tmp), colData = tmp, design = ~ comments)
sfs <- log10(colSums(counts[, colnames(counts_tmp)])) / median(log10(colSums(counts[, colnames(counts_tmp)])))
sizeFactors(des) <- sfs
sizeFactors(des) <- rep(1, ncol(des))
des <- DESeq(des)
res <- results(des, contrast = c("comments", "touched by users", "untouched"))
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res$pvalue),]

ggplot(res, aes(log2FoldChange, -log10(pvalue), color = Order)) +
  geom_point() +
  theme_bw()

# Create some plots ####
tmp <- log(OTU@.Data + 1)
tmp <- tmp[which(apply(tmp, 1, var) > 0),]
pca <- prcomp(t(tmp), scale. = T)

counts <- OTU@.Data

aframe <- data.frame(pca$x, data.frame(info), total_counts = colSums(counts))
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

