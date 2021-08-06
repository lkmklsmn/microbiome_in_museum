# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)

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

ok = is.na(as.character(TAX[,"Phylum"]))
GP1 = prune_taxa(ok, p)

ok = names(which(colSums(OTU) > 100))
GP1 = prune_samples(ok, p)

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(physeq = GP1, method = "MDS")
plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
plot_ordination(GP1, GP.ord, type = "samples", color = "Object.number", shape = 'Institution', title = "samples") 

# Plot PCA ####
aframe <- data.frame(GP.ord$vectors, sample_data(GP1))
aframe$label <- NA
aframe$label[aframe$type.of.room == ""] <- aframe$sample.id[aframe$type.of.room == ""]

ggplot(aframe, aes(Axis.1, Axis.2, shape = institution, color = object, label = label)) +
  geom_point() + ggrepel::geom_label_repel() +
  theme_bw()

ggplot(aframe, aes(Axis.1, Axis.2, shape = institution, color = spot == "Air")) +
  geom_point() +
  theme_bw()

ggplot(aframe, aes(Axis.3, Axis.4, shape = institution, color = factor(object_number))) +
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

ggplot(aframe, aes(X1, X2, color = institution, shape = institution, label = sample.id)) +
  geom_point() + geom_label() +
  theme_classic()

ggplot(aframe, aes(X1, X2, color = log10(total), shape = institution)) +
  geom_point() +
  theme_classic()

ggplot(aframe, aes(PC1, PC2, color = log10(total), shape = institution)) +
  geom_point() +
  theme_classic()

ggplot(aframe, aes(PC1, PC2, color = institution, shape = institution, label = sample.id)) +
  geom_point() + geom_label() +
  theme_classic()

# Run positive unlabeled learning ####
outcome <- meta$institution == "DSMZ"
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

p1 <- ggplot(aframe, aes(avg_pred, log10(total), color = collection == "Control", shape = institution)) +
  geom_point() +
  ylab("Total count [log10]") + xlab("Average prediction") +
  ggtitle("Positive-unlabeled learning") +
  theme_bw()

p2 <- ggplot(aframe, aes(avg_pred, log10(total), color = spot == "Air", shape = institution)) +
  geom_point() +
  ylab("Total count [log10]") + xlab("Average prediction") +
  ggtitle("Positive-unlabeled learning") +
  theme_bw()

gridExtra::grid.arrange(p1, p2, ncol = 2)

# Run analysis on good lion samples ####
ok <- which(aframe$avg_pred < 1.5)
tmp <- meta[ok,]
ok_samples <- tmp$sample.id[grep("Lion", tmp$object)]
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
aframe_tmp <- data.frame(meta[match(ok_samples, meta$sample.id), ], total = colSums(counts_tmp), mds)

ggplot(aframe_tmp, aes(X1, X2, color = comments, shape = object)) +
  geom_point() +
  ggtitle("MDS - Lion") +
  theme_classic()

# Hellenistic Hall ####
ok_samples <- tmp$sample.id[grep("Hellenistic", tmp$object)]
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
tmp <- meta[which(meta$object == "Tendaguru"), ]
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
  xlab("MDS1") + ylab("MDS2") +
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

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Streptococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Streptococcus abundance")

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Staphylococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Streptococcus abundance")

# Mollucs ####
tmp <- meta[which(meta$collection == "Mollucs"), ]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = comments, shape = object)) +
  geom_point() +
  ggtitle("MDS - Mollucs") +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = comments, shape = object)) +
  geom_point() +
  ggtitle("PCA - Mollucs") +
  theme_classic()

aframe_tmp <- data.frame(PC1 = pca$rotation[,1], TAX@.Data[rownames(pca$rotation), ])
asplit <- split(aframe_tmp$PC1, aframe_tmp$Genus)
aframe_tmp$Genus <- factor(aframe_tmp$Genus, levels = names(sort(unlist(lapply(asplit, median)))))
ggplot(aframe_tmp[!is.na(aframe_tmp$Genus),], aes(PC1, Genus, color = Phylum)) +
  geom_boxplot() +
  xlab("PC1 loading") +
  theme_bw()

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Streptococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Streptococcus abundance")

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Staphylococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Staphylococcus abundance")


# Combine into one analysis ####
tmp <- meta[meta$institution == "MfN",]
tmp <- tmp[tmp$comments != "reference",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

aframe_tmp$touched <- "yes"
aframe_tmp$touched[grep("untouched", aframe_tmp$comments)] <- "no"
ggplot(aframe_tmp, aes(X1, X2, color = touched, shape = collection)) +
  geom_point() +
  ggtitle("MDS - Mollucs") +
  xlab("MDS1") + ylab("MDS2") +
  theme_classic()
ggsave("")

ggplot(aframe_tmp, aes(PC1, PC2, color = touched, shape = object)) +
  geom_point() +
  ggtitle("MDS - Mollucs") +
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

anno_col <- data.frame(aframe_tmp[, c("touched", "collection", "object", "avg_pred")])
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(TAX@.Data[rownames(normalized),])
rownames(anno_row) <- rownames(normalized)

sig <- rownames(res)[p.adjust(res[,2], method = "BH") < 0.1]
tmp <- normalized[sig, anno_col$collection == "Mollucs"]

pheatmap(normalized[sig, order(anno_col$collection, anno_col$touched, anno_col$object)],
         #scale = "row",
         cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = 10,
         annotation_col = anno_col,
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")])

plot_feature <- function(feature){
  aframe <- data.frame(expr = normalized[feature, ], aframe_tmp[, c("touched", "object", "collection")])
  nom <- TAX@.Data[feature, -1]
  nom <- paste(nom, collapse = ",")
  ggplot(aframe, aes(touched, expr, color = touched)) +
    facet_wrap(~ collection) +
    geom_boxplot() + geom_point() +
    ggtitle(nom) +
    ylab("Abundance") +
    theme_bw()
}
plot_feature("0664a0c96e567c98809225119f8fad99")

tmp <- data.frame(feature = colSums(counts[sig,]), meta)
tmp$touched <- "no"
tmp$touched[tmp$comments == "touched"] <- "yes"
ggplot(tmp[grep("Lion", tmp$object),], aes(touched, log10(feature + 1), color = touched)) +
  geom_boxplot() + geom_point() +
  ggtitle("Gate Lion") +
  ylab("'Touch' signature") +
  theme_bw()

# SPK ####
tmp <- meta[meta$institution == "SPK",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = type_of_room, shape = type_of_room)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  xlab("MDS1") + ylab("MDS2") +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = type_of_room, shape = material)) +
  geom_point() +
  ggtitle("PCA - SPK") +
  theme_classic()

# SPK Lions ####
tmp <- meta[grep("Lion", meta$object),]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = comments, shape = object)) +
  geom_point() +
  ggtitle("MDS - Gate Lion") +
  xlab("MDS1") + ylab("MDS2") +
  theme_classic()

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Streptococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Streptococcus abundance")

ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == "g__Staphylococcus")))
boxplot(split(colMeans(normalized[ok_features,]), tmp$comments),
        ylab = "Staphylococcus abundance")


# ####
tmp <- meta[meta$type_of_room == "exhibition",]
tmp <- tmp[tmp$object != "",]
tmp <- tmp[tmp$object %in% c("Ishtar Gate, right side", "Procession Street, left side") ,]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = object, shape = object)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()

ggplot(aframe_tmp, aes(X1, X2, color = avg_pred, shape = material)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()


ggplot(aframe_tmp, aes(PC1, PC2, color = object, shape = object)) +
  geom_point() +
  ggtitle("PCA - SPK") +
  theme_classic()

aframe_tmp$condition <- "yes"
aframe_tmp$condition[grep("white", aframe_tmp$spot)] <- "no"
ggplot(aframe_tmp, aes(PC1, PC2, color = condition, shape = object)) +
  geom_point() +
  ggtitle("PCA - SPK") +
  theme_classic()

# Procession Street ####
tmp <- meta[meta$object == "Procession Street, left side",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

aframe_tmp$height <- c("low", "low", "middle", "middle", "high", "high", "very high", "very high")
ggplot(aframe_tmp, aes(X1, X2, color = height, shape = object)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()

ggplot(aframe_tmp, aes(X1, X2, color = avg_pred, shape = object)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()


# Procession Street ####
tmp <- meta[meta$object == "Ishtar Gate, right side",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = spot, shape = object)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = spot, shape = object)) +
  geom_point() +
  ggtitle("PCA - SPK") +
  theme_classic()


# Zeus Sosipolis Tempel, Hellenistic Hall ####
tmp <- meta[meta$object == "Zeus Sosipolis Tempel, Hellenistic Hall",]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[which(zeros > 1),]

counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

mds <- cmdscale(dist(t(normalized)), k = 2)
pca <- prcomp(t(normalized), scale. = T)
aframe_tmp <- data.frame(meta[match(colnames(counts_tmp), meta$sample.id), ], total = colSums(counts_tmp), mds, pca$x)

ggplot(aframe_tmp, aes(X1, X2, color = spot, shape = object)) +
  geom_point() +
  ggtitle("MDS - SPK") +
  theme_classic()

# Create some plots ####
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
