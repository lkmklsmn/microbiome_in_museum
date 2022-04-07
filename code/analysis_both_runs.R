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

# Load sample sheet ####
meta <- readxl::read_excel("/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/sample_sheet_data_science.xlsx", sheet = "Tabelle1")
meta <- data.frame(meta)
rownames(meta) <- meta$sample.id

# Load taxonomy info ####
classified <- read_qza("/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/taxonomy.qza")
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

# Plot total counts ####
counts <- read_qza("/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/table.qza")$data
total <- colSums(counts)
num_detected <- apply(counts, 2, function(x) sum(x > 0))

meta <- data.frame(meta[match(colnames(counts), meta[[1]]), ])

aframe <- data.frame(total, num_detected, meta)
ggplot(aframe, aes(total + 1, num_detected,
                   label = sample.id,
                   color = interaction(Sequencing.run, Institution))) +
  geom_point() + geom_label() +
  scale_x_continuous(trans = "log10") +
  ylab("Number of feature detected") +
  xlab("Total number of reads") +
  theme_bw()

zeros <- apply(counts, 1, function(x) sum(x > 0))

tmp <- counts[which(zeros > 1),]
normalized <- t(t(tmp)/colSums(tmp))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe <- data.frame(meta, total = colSums(counts), pca$x[, 1:10])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

ggplot(aframe, aes(PC1, PC2, color = Institution, shape = Institution, label = sample.id)) +
  geom_point() + geom_label() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()

# Create phylum barplot ####
level <- "Phylum"
ok <- which(colSums(counts) > 200)
tmp <- t(t(counts[, ok])/colSums(counts[, ok]))
asplit <- split(rownames(TAX), TAX[, level])
asplit <- asplit[which(unlist(lapply(asplit, length)) > 1)]
norm <- do.call(rbind, lapply(asplit, function(x) colSums(counts[x, ok])))
norm <- t(t(norm)/colSums(norm))

if(level=="Genus") top_hits <- tail(setdiff(names(sort(rowMeans(norm))), "g__"), 15)
if(level=="Phylum") top_hits <- tail(setdiff(names(sort(rowMeans(norm))), "p__"), 10)
rest <- colSums(norm[setdiff(rownames(norm), top_hits),])

aframe <- data.frame(t(norm[top_hits, ]), rest, meta[ok, ])
ord <- order(aframe$Institution, aframe$Object)
aframe$sample.id <- factor(aframe$sample.id, levels = (aframe$sample.id[ord]))

aframe <- reshape2::melt(aframe, id.vars = colnames(meta))

farben <- qualitative_hcl(length(top_hits) + 1, "Dark2")
farben[length(farben)] <- "grey"

ggplot(aframe,aes(x = sample.id, y = value, fill = variable)) + 
  facet_grid(Sequencing.run ~ Institution, scales = "free_x") +
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_manual(values = farben, name = "Phylum") +
  ylab("Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Combine into one analysis ####
tmp <- meta[which(meta$Institution == "MfN" & meta$Comments != "reference"),]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros > 2)

tmp <- meta[which(meta$Institution == "MfN" & meta$Sequencing.run == "first" & meta$Comments != "reference"),]
ok_samples_1st <- tmp$sample.id
rownames(tmp) <- ok_samples_1st
counts_tmp <- counts[rows_ok, ok_samples_1st]
#zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
counts_tmp <- counts_tmp[, colSums(counts_tmp) > 100]
ok_samples_1st <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

rows_ok <- names(which(apply(normalized, 1, var) > 0))
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

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing.run == "second" & meta$Comments != "reference"),]
ok_samples_2nd <- tmp$sample.id
counts_tmp <- counts[, ok_samples_2nd]
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 200]
ok_samples_2nd <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

preds <- predict(pca, t(normalized))
preds_seq2 <- predict(lda, t(normalized))

aframe_tmp <- data.frame(rbind(pca$x[, 1:5], preds[, 1:5]),
                         meta[c(rownames(pca$x), rownames(preds)), ])

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing.run)) +
  geom_point() +
  labs(title = "PCA - Tendaguru", subtitle = "2nd sequencing run does not separate touched from untouched") +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  scale_color_manual(values = c("red", "lightgrey", "grey","darkgrey", "darkred")) +
  theme_classic()

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing.run)) +
  facet_wrap(~ Sequencing.run, scales ="free") +
  geom_point() +
  labs(title = "PCA - Tendaguru", subtitle = "2nd sequencing run does not separate touched from untouched") +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  scale_color_manual(values = c("red", "lightgrey", "grey","darkgrey", "darkred")) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic()

boxplot(split(preds_seq1$x, meta[rownames(preds_seq1$x), "Comments"]), las = 2)
boxplot(split(preds_seq2$x, meta[rownames(preds_seq2$x), "Comments"]))

# Tendaguru runs 1 & 2 ####
tmp <- meta[which(meta$Object == "Tendaguru" & meta$Comments != "reference"),]
ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros > 1)

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing.run == "first" & meta$Comments != "reference"),]
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

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing.run == "second" & meta$Comments != "reference"),]
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

ggplot(aframe_tmp, aes(PC1, PC2, color = Comments, shape = Sequencing.run)) +
  facet_wrap(~ Sequencing.run, scales = "free") + 
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

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing.run == "first" & meta$Comments != "reference"),]
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

tmp <- meta[which(meta$Object == "Tendaguru" & meta$Sequencing.run == "second" & meta$Comments != "reference"),]
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
tmp <- meta[which(meta$Institution == "MfN" & meta$Sequencing.run == "first" & meta$Comments != "reference"),]
ok_samples_1st <- tmp$sample.id
rownames(tmp) <- ok_samples_1st
counts_tmp <- counts[, ok_samples_1st]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(apply(counts_tmp, 1, function(x) sum(x > 5)) > 2)
counts_tmp <- counts_tmp[rows_ok, colSums(counts_tmp) > 100]
ok_samples_1st <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)
rows_ok <- names(which(apply(normalized, 1, var) > 0))
normalized <- normalized[rows_ok,]

res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x, meta[colnames(normalized), c("Comments", "Object", "Collection")])
  aframe$touched <-"yes"
  aframe$touched[grep("untouched", aframe$Comments)] <- "no"
  afit <- try(summary(aov(expr ~ Collection + touched, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  as.numeric(afit[[1]][2, 4:5])
}))
colnames(res) <- c("F_value", "p_value")
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res[, 2]),]
res$p_adjusted <- p.adjust(res$p_value, method = "BH")

anno_col <- meta[colnames(normalized), c("Comments", "Collection", "Object")]
anno_col$touched <-"yes"
anno_col$touched[grep("untouched", anno_col$Comments)] <- "no"
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(TAX@.Data[rownames(normalized),])
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
                  meta,
                  total)

tmp$touched <- "yes"
tmp$touched[grep("untouched", tmp$Comments)] <- "no"
tmp <- tmp[grep("Lion", tmp$Object),]
tmp <- tmp[tmp$Sequencing.run == "first",]

summary(lm(feature_norm ~ total + touched, data = tmp[tmp$Sequencing.run == "first",]))

ggplot(tmp, aes(touched, log10(feature + 1), color = touched)) +
  facet_wrap(~Object) +
  geom_boxplot() + geom_point() +
  ggtitle("Gate Lion") +
  ylab("'Touch' signature") +
  theme_bw()

plot_genus <- function(genus){
  ok_features <- intersect(rownames(normalized), names(which(TAX@.Data[, "Genus"] == genus)))
  
  aframe <- data.frame(tmp,
                       feature = colMeans(normalized[ok_features, rownames(tmp)]))
  
  ggplot(aframe, aes(touched, feature, color = touched)) +
    facet_wrap(~Object) +
    geom_boxplot() + geom_point() +
    ggtitle("Gate Lion") +
    ylab(genus) +
    theme_bw()
  
}
plot_genus("g__Streptococcus")
plot_genus("g__Staphylococcus")

# Define functions ####
get_subset <- function(samples){
  ok_samples_2nd <- samples
  
  counts_tmp <- counts[, ok_samples_2nd]
  counts_tmp <- counts_tmp[which(rowSums(counts_tmp) > 10), which(colSums(counts_tmp) > 100)]
  
  ok_samples_2nd <- colnames(counts_tmp)
  
  normalized <- t(t(counts_tmp)/colSums(counts_tmp))
  normalized <- sqrt(normalized)
  
  pca <- prcomp(t(normalized), scale. = T)
  var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)
  
  aframe_tmp <- data.frame(pca$x[, 1:5],
                           meta[ok_samples_2nd, ])
  aframe_tmp
}

run_pca <- function(samples, color = "Object", label = "sample.id"){
  
  aframe_tmp <- get_subset(samples = samples)
  colnames(aframe_tmp) <- gsub(color, "color", colnames(aframe_tmp))
  colnames(aframe_tmp) <- gsub(label, "label", colnames(aframe_tmp))
  
  ggplot(aframe_tmp, aes(PC1, PC2, color = color, label = label)) +
    geom_point() +
    ggrepel::geom_label_repel() +
    ggtitle("PCA - SPK") +
    xlab(paste("PC1", var_explained[1], "%")) +
    ylab(paste("PC2", var_explained[2], "%")) +
    theme_classic()
}

normalize_counts <- function(counts){
  counts <- counts[rowSums(counts) > 0,]
  sqrt(t(t(counts)/colSums(counts))) 
}

# SPK analysis - height profile ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Collection == "VAM" & 
                   meta$Sequencing.run == "second" &
                   meta$Material == "Glazed Ceramic" &
                   is.na(meta$Comments)), "sample.id"]
run_pca(samples = ok, color = "Height", label = "Object")

aframe_tmp <- get_subset(samples = ok)
ggplot(aframe_tmp, aes(PC1, as.numeric(Height), color = Height, label = sample.id)) +
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

# SPK analysis -  Glazed Ceramic vs Stone ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Type.of.room == "exhibition" &
                   is.na(meta$Comments)), "sample.id"]
run_pca(samples = ok, color = "Material", label = "Object")

# SPK analysis - werkstatt vs exhibition ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Material == "Glazed Ceramic"), "sample.id"]
run_pca(samples = ok, color = "Type.of.room", label = "Object")

# SPK analysis - blue vs yellow tiles ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Type.of.room == "exhibition" &
                   meta$Material == "Glazed Ceramic" &
                   is.na(meta$Comments)), "sample.id"]
run_pca(samples = ok, color = "Color", label = "Object")

# SPK analysis - STS vs RRP ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Comments %in% c("StS", "visitors and StS", "RRP")), "sample.id"]
run_pca(samples = ok, color = "Comments", label = "Object")
run_pca(samples = ok, color = "Material", label = "Object")

# Touched vs untouched ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
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


