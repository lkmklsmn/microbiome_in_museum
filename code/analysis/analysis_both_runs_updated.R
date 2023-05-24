# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(microbiome)
library(colorspace)

#setwd("/Users/lukas/OneDrive/Documents/GitHub/lounsbery/")
# Load data ####
load("../../data/Touched_DADA2_objects_2022-07-28.RData")
counts <- data.matrix(seqtab.final[, rownames(metadata)])
rownames(counts) <- seqtab.final$Sequence

# Plot total counts ####
tmp <- metadata
lvs <- tmp$sample_id[order(tmp$Sequencing_run, tmp$Institution, tmp$Collection, -tmp$Seq_depth)]
tmp$sample_id <- factor(tmp$sample_id, levels = lvs)
ggplot(tmp, aes(sample_id, Seq_depth, color = Sequencing_run)) +
  facet_wrap(~ Sequencing_run ~ Institution ~ Collection, scales = "free_x", nrow = 1) +
  geom_point() + 
  scale_y_continuous(trans = "log10") +
  theme_bw() + theme(axis.text.x = element_blank())

tmp$detected <- apply(counts, 2, function(x) sum(x > 0))
ggplot(tmp, aes(Seq_depth, detected, color = Sequencing_run, shape = Institution)) +
  geom_point() + 
  labs(x = "Total reads",
       y = "Number of ASVs detected") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_bw()

# Plot PCA of all samples ####
tmp <- counts[, colSums(counts) > 0]
tmp <- tmp[which(apply(tmp, 1, var) > 0), ]
normalized <- t(t(tmp)/colSums(tmp))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe <- data.frame(metadata[colnames(normalized), ],
                     pca$x[, 1:10])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

ggplot(aframe, aes(PC1, PC2,
                   color = Institution,
                   label = sample_id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

ggplot(aframe, aes(PC1, PC2,
                   color = Sequencing_run,
                   label = sample_id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

# Combine into one analysis ####
tmp <- metadata[which(metadata$Institution == "MfN" &
                        metadata$Comments != "reference"),]
ok_samples <- tmp$sample_id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros > 2)

tmp <- metadata[which(metadata$Institution == "MfN" &
                        metadata$Sequencing_run == "first" &
                        metadata$Comments != "reference"),]
ok_samples_1st <- tmp$sample_id
counts_tmp <- counts[rows_ok, ok_samples_1st]
#zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

rows_ok <- which(apply(normalized, 1, var) > 0)
normalized <- normalized[rows_ok,]

lda <- MASS::lda(t(normalized), grouping = tmp[colnames(normalized), "Comments"] != "untouched")
preds_seq1 <- predict(lda, t(normalized))

pca <- prcomp(t(normalized))
aframe <- data.frame(pca$x, tmp[colnames(normalized), ])
aframe$touched <- "yes"
aframe$touched[grep("untouched", aframe$Comments)] <- "no"
ggplot(aframe, aes(PC1, PC2, color = touched, shape = Collection)) +
  labs(title = "MfN: Confirmatory results", subtitle = "Touched molluscs and fossils cluster together") +
  geom_point() +
  theme_classic()

tmp <- metadata[which(metadata$Object == "Tendaguru" &
                        metadata$Sequencing_run == "second" &
                        metadata$Comments != "reference"),]
counts_tmp <- counts[, tmp$sample_id]
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 200]
ok_samples_2nd <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

preds <- predict(pca, t(normalized))
preds_seq2 <- predict(lda, t(normalized))

aframe_tmp <- data.frame(rbind(pca$x[, 1:5], preds[, 1:5]),
                         metadata[c(rownames(pca$x), rownames(preds)), ])

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing_run)) +
  geom_point() +
  labs(title = "PCA - Tendaguru", subtitle = "2nd sequencing run does not separate touched from untouched") +
  scale_color_manual(values = c("red", "lightgrey", "grey","darkgrey", "darkred")) +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing_run)) +
  facet_wrap(~ Sequencing_run, scales ="free") +
  geom_point() +
  labs(title = "PCA - Tendaguru",
       subtitle = "2nd sequencing run does not separate touched from untouched") +
  scale_color_manual(values = c("red", "lightgrey", "grey","darkgrey", "darkred")) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic()

lvs <- c("very much touched by users",
         "touched by users",
         "untouched",
         "untouched, original tin unopened",
         "untouched, original package unopened")
aframe_tmp$Comments <- factor(aframe_tmp$Comments, levels = lvs)
ggplot(aframe_tmp, aes(Comments, PC1, fill = Comments)) +
  facet_wrap(~ Sequencing_run, scales = "free_x") +
  labs(title = "PCA - Tendaguru",
       subtitle = "PC1 separates touched vs untouched objects in 1st but not 2nd sequencing run") +
  geom_boxplot() +
  scale_fill_manual(values = c("darkred", "red", "lightgrey", "grey","darkgrey")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Tendaguru runs 1 & 2 ####
tmp <- meta[which(meta$Object == "Tendaguru" & meta$Comments != "reference"),]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros > 1)

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing_run == "first" & meta$Comments != "reference"),]
ok_samples_1st <- tmp$sample.id
counts_tmp <- counts[, ok_samples_1st]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_1st <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe <- data.frame(pca$x, tmp)
ggplot(aframe, aes(PC1, PC2, color = Inv..number)) +
  geom_point()

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing_run == "second" & meta$Comments != "reference"),]
ok_samples_2nd <- tmp$sample.id
counts_tmp <- counts[, ok_samples_2nd]
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_2nd <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

preds <- predict(pca, t(normalized))

aframe_tmp <- data.frame(rbind(pca$x[, 1:5], preds[, 1:5]),
                         meta[match(c(ok_samples_1st, ok_samples_2nd), meta$sample.id), ])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing_run)) +
  facet_wrap(~ Sequencing_run, scales = "free") + 
  geom_point() +
  labs(title = "PCA - Tendaguru", subtitle = "2nd sequencing run does not separate touched from untouched") +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  theme_classic()

# Run LDA model ####
tmp <- meta[which(meta$Object == "Tendaguru" & meta$Comments != "reference"),]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros > 1)

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing_run == "first" & meta$Comments != "reference"),]
ok_samples_1st <- tmp$sample.id
rownames(tmp) <- ok_samples_1st
counts_tmp <- counts[, ok_samples_1st]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_1st <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)
rows_ok <- names(which(apply(normalized, 1, var) > 0))
normalized <- normalized[rows_ok,]

lda <- MASS::lda(t(normalized), grouping = tmp[colnames(normalized), "Comments"])
preds_seq1 <- predict(lda, t(normalized))

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing_run == "second" & meta$Comments != "reference"),]
ok_samples_2nd <- tmp$sample.id
rownames(tmp) <- ok_samples_2nd
counts_tmp <- counts[, ok_samples_2nd]
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_2nd <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

preds_seq2 <- predict(lda, t(normalized))
boxplot(split(preds_seq1$x, meta[rownames(preds_seq1$x), "Comments"]), ylab = "LD1", main = "Seq run 1")
boxplot(split(preds_seq2$x, meta[rownames(preds_seq2$x), "Comments"]), ylab = "LD1", main = "Seq run 2")

# Run DE ####
tmp <- metadata[which(metadata$Institution == "MfN" &
                        metadata$Sequencing_run == "first" &
                        metadata$Comments != "reference"),]
counts_tmp <- counts[, rownames(tmp)]

zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(apply(counts_tmp, 1, function(x) sum(x > 5)) > 2)
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_1st <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)
rows_ok <- which(apply(normalized, 1, var) > 0)
normalized <- normalized[rows_ok,]

res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x,
                       metadata[colnames(normalized),
                                c("Comments", "Object", "Collection")])
  aframe$touched <-"yes"
  aframe$touched[grep("untouched", aframe$Comments)] <- "no"
  afit <- try(summary(aov(expr ~ Collection + touched, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  as.numeric(afit[[1]][2, 4:5])
}))
colnames(res) <- c("F_value", "p_value")
cols <- c("Kingdom", "Phylum", "Class", "Order", "Family",
          "Genus", "IdTaxa_RDP", "IdTaxa_Silva", "Mixed_genus")
res <- data.frame(res, seqtab.final[match(rownames(res), seqtab.final$Sequence), cols])
res <- res[sort.list(res[, 2]),]
res$p_adjusted <- p.adjust(res$p_value, method = "BH")

anno_col <- metadata[colnames(normalized), c("Comments", "Collection", "Object")]
anno_col$touched <-"yes"
anno_col$touched[grep("untouched", anno_col$Comments)] <- "no"
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(seqtab.final[match(rownames(normalized), seqtab.final$Sequence), cols])
rownames(anno_row) <- rownames(normalized)

sig <- rownames(res)[res$p_adjusted < 0.1]

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollucs` = "blue", `VertPal` = "darkblue")
object <- c(`Chiton` = "bisque3", `Tendaguru` = "darksalmon", `Triton horn`= "yellow", `Mobbi collection`= "brown")
ann_colors <- list(touched = touched,
                   collection = collection,
                   object = object)

pheatmap(normalized[sig, order(anno_col$Collection, anno_col$touched, anno_col$Object)],
         cluster_cols = F,
         gaps_col = 11,
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col[, c("Object", "Collection", "touched")],
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")],
         annotation_colors = ann_colors)

sig_touch_1st <- sig

# Validation of touch signature in Gate lions ####
normalized <- t(t(counts)/colSums(counts))
normalized <- sqrt(normalized)

tmp <- data.frame(feature = colSums(counts[sig,]),
                  feature_norm = colMeans(normalized[sig,]),
                  detected = apply(counts[sig, ], 2, function(x) sum(x > 0)),
                  metadata)

tmp$touched <- "yes"
tmp$touched[grep("untouched", tmp$Comments)] <- "no"
tmp <- tmp[grep("Lion", tmp$Object),]
tmp <- tmp[tmp$Sequencing_run == "first",]

summary(glm(feature ~ log(Seq_depth) + touched, data = tmp, family = poisson))

summary(MASS::glm.nb(feature ~ log(Seq_depth) + touched, data = tmp))

ggplot(tmp, aes(touched, log10(feature + 1), color = touched)) +
  facet_wrap(~Object) +
  geom_boxplot() + geom_point() +
  labs(title = "Gate Lion",
       subtitle = "Gate Lions confirm 'touch' signature") +
  ylab("'Touch' signature") +
  theme_bw()

plot_genus <- function(genus){
  ok_features <- intersect(rownames(normalized),
                           seqtab.final$Sequence[which(seqtab.final$Genus == genus)])
  
  aframe <- data.frame(tmp,
                       feature = colMeans(normalized[ok_features,
                                                     rownames(tmp)]))
  
  ggplot(aframe, aes(touched, feature, color = touched)) +
    facet_wrap(~ Object) +
    geom_boxplot() + geom_point() +
    ggtitle("Gate Lion") +
    ylab(genus) +
    theme_bw()
}
plot_genus("Streptococcus")
plot_genus("Staphylococcus")

# Define functions ####
get_subset <- function(samples){
  
  counts_tmp <- counts[, samples]
  counts_tmp <- counts_tmp[which(rowSums(counts_tmp) > 10), which(colSums(counts_tmp) > 100)]
  
  ok_samples_2nd <- colnames(counts_tmp)
  
  normalized <- t(t(counts_tmp)/colSums(counts_tmp))
  normalized <- sqrt(normalized)
  
  pca <- prcomp(t(normalized), scale. = T)
  var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)
  
  aframe_tmp <- data.frame(pca$x[, 1:5],
                           metadata[ok_samples_2nd, ])
  aframe_tmp
}

run_pca <- function(samples, color = "Object", label = "sample.id"){
  
  aframe_tmp <- get_subset(samples = samples)
  colnames(aframe_tmp) <- gsub(color, "color", colnames(aframe_tmp))
  colnames(aframe_tmp) <- gsub(label, "label", colnames(aframe_tmp))
  
  ggplot(aframe_tmp, aes(PC1, PC2, color = color, label = label)) +
    geom_point() +
    ggrepel::geom_label_repel() +
    labs(title = "PCA - SPK") +
    theme_classic()
}

normalize_counts <- function(counts){
  counts <- counts[rowSums(counts) > 0,]
  sqrt(t(t(counts)/colSums(counts))) 
}

# SPK analysis - height profile ####
ok <- metadata[which(metadata$Institution == "SPK" &
                       metadata$Collection == "VAM" & 
                       metadata$Sequencing_run == "second" &
                       metadata$Material == "Glazed Ceramic" &
                   is.na(metadata$Comments)), "sample_id"]
run_pca(samples = ok, color = "Height", label = "Object")

aframe_tmp <- get_subset(samples = ok)
ggplot(aframe_tmp, aes(PC1, as.numeric(Height), color = Height, label = sample_id)) +
  geom_point() + geom_label() +
  ylab("Height [m]") +
  xlab("PC1") +
  theme_classic()
summary(lm(PC1 ~ as.numeric(Height), data = aframe_tmp))

# Create heatmap of ASVs associated with height #####
normalized <- normalize_counts(counts[,ok])
res_height <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x,
                       height = as.numeric(meta[colnames(normalized), c("Height")]))
  
  afit <- try(summary(lm(expr ~ height, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  coefficients(afit)[2, c(1,4)]
}))
res_height <- data.frame(res_height)
colnames(res_height) <- c('estimate', 'pval')

res_height <- res_height[sort.list(res_height$pval),]
res_height$p_adjusted <- p.adjust(res_height$pval, method = "BH")

anno_col <- meta[colnames(normalized), ]
anno_col$Height <- as.numeric(anno_col$Height)
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(TAX@.Data[rownames(normalized),])
rownames(anno_row) <- rownames(normalized)

sig <- rownames(res_height)[1:20]

pheatmap(normalized[sig, order(anno_col$Height)],
         cluster_cols = F,
         scale = "row", breaks = seq(-1.5, 1.5, length = 100),
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col[, c("Height", "Object", "Spot")],
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")])


# SPK analysis - Gate lions (both runs) ####
ok <- metadata[which(metadata$Institution == "SPK" &
                       metadata$Object %in% c("Gate Lion left Sam'al", "Gate Lion right Sam'al") &
                       metadata$Sequencing_run == "first"),
                "sample_id"]
aframe_tmp <- get_subset(samples = ok)

ggplot(aframe_tmp, aes(PC1, PC2, color = Object, label = Comments)) +
  facet_wrap(~ Sequencing_run, scales = "free", nrow = 2) +
  geom_point() + geom_label() +
  theme_classic()

# SPK analysis -  Glazed Ceramic vs Stone ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing_run == "second" &
                   meta$Type.of.room == "exhibition" &
                   is.na(meta$Comments)), "sample.id"]
run_pca(samples = ok, color = "Material", label = "Object")

# SPK analysis - werkstatt vs exhibition ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing_run == "second" &
                   meta$Material == "Glazed Ceramic"), "sample.id"]
run_pca(samples = ok, color = "Type.of.room", label = "Object")

# SPK analysis - blue vs yellow tiles ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing_run == "second" &
                   meta$Type.of.room == "exhibition" &
                   meta$Material == "Glazed Ceramic" &
                   is.na(meta$Comments)), "sample.id"]
run_pca(samples = ok, color = "Color", label = "Object")

# SPK analysis - STS vs RRP ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing_run == "second" &
                   meta$Comments %in% c("StS", "visitors and StS", "RRP")), "sample.id"]
run_pca(samples = ok, color = "Comments", label = "Object")
run_pca(samples = ok, color = "Material", label = "Object")

# Touched vs untouched ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing_run == "second" &
                   meta$Type.of.room == "exhibition" &
                   meta$sample.id != "2Pergamon Mus. Negativprobe"), "sample.id"]
run_pca(samples = ok, color = "Comments", label = "Object")

aframe_tmp <- get_subset(samples = ok)
aframe_tmp$touch <- "no"
aframe_tmp$touch[!is.na(aframe_tmp$Comments)] <- "yes"
ggplot(aframe_tmp, aes(Comments, PC1, color = touch, label = sample.id)) +
  geom_boxplot() + geom_point() +
  labs(title = "SPK", subtitle = "Touched vs untouched", y = "PC1", x = "Comments") +
  theme_classic()
summary(lm(PC1 ~ touch, data = aframe_tmp))

# Validate touch signature from 1st run ####
normalized <- normalize_counts(counts[, ok])

means <- colMeans(normalized[intersect(rownames(normalized), sig_touch_1st),])
aframe_tmp$mean <- means

sums <- colSums(counts[sig_touch_1st, ok])
aframe_tmp$sum <- sums

ggplot(aframe_tmp, aes(Comments, mean, color = touch, label = sample.id)) +
  geom_boxplot() + geom_point() +
  labs(title = "SPK", subtitle = "Touched vs untouched",
       y = "Touch signature (Tendaguru)", x = "Comments") +
  theme_classic()

ggplot(aframe_tmp, aes(Comments, sum, color = touch, label = sample.id)) +
  geom_boxplot() + geom_point() +
  labs(title = "SPK", subtitle = "Touched vs untouched",
       y = "Touch signature (Tendaguru)", x = "Comments") +
  theme_classic()

# ####


