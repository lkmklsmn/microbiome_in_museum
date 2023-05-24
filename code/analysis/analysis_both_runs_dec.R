# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
#library(microbiome)
#library(colorspace)

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

# Load ASV count table ####
counts <- read_qza("/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/table.qza")$data
counts <- counts[which(apply(counts, 1, function(x) sum(x > 0)) >= 2), ]
tax_info <- TAX@.Data[rownames(counts), ] 

# Calculate total counts & number features detected ####
meta <- data.frame(meta[match(colnames(counts), meta$sample.id), ])
meta$total_reads <- colSums(counts)
meta$num_detected <- apply(counts, 2, function(x) sum(x > 0))

# Plot total counts ####
aframe <- data.frame(meta)
aframe$Institution[aframe$Collection == "RF"] <- "RF"
aframe$group <- paste(aframe$Institution)
aframe$group[aframe$group == "DSMZ"] <- "Controls"

aframe$group[aframe$Sampling.number %in% c("Pergamonmuseum  Negativprobe",
                                           "RF Negativprobe",
                                           "2VAM NC Werkstatt")] <- "Controls"

ggplot(aframe, aes(total_reads + 1, num_detected,
                   color = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 500, linetype = "dashed", color = "grey") +
  geom_point() + 
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c("darkgrey", "darkgreen", "darkblue", "darkred")) +
  labs(title = "", subtitle = "",
       y = "Number of feature detected", x = "Total number of reads") +
  theme_classic()
ggsave("../../figs/total_counts_num_detected_all_samples.pdf",
       height = 6, width = 7)

# Remove sample w low counts or ASV features detected ####
bad <- which(aframe$num_detected < 100 | aframe$total_reads < 500 | aframe$group == "Controls")
meta <- meta[-bad, ]
counts <- counts[, -bad]

# Plot PCA ####
zeros <- apply(counts, 1, function(x) sum(x > 0))

tmp <- counts[which(zeros > 1),]
normalized <- t(t(tmp)/colSums(tmp))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
aframe <- data.frame(meta,
                     total = colSums(counts), pca$x[, 1:10])

var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)

ggplot(aframe, aes(PC1, PC2,
                   color = Sequencing.run,
                   shape = Institution)) +
  geom_point() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  ggtitle("Principal component analysis") +
  theme_classic()
ggsave("../../figs/pca_all_samples_after_filtering.pdf",
       height = 6,
       width = 7)

# Create phylum barplot ####
level <- "Phylum"
ok <- which(colSums(counts) > 200)
tmp <- t(t(counts[, ok])/colSums(counts[, ok]))
asplit <- split(rownames(tax_info), tax_info[, level])
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
aframe$variable <- gsub("p__", "", fixed = T, aframe$variable)

farben <- colorspace::qualitative_hcl(length(top_hits) + 1, "Dark2")
farben[length(farben)] <- "grey"

aframe$label <- paste(aframe$Object, aframe$Sampling.number)
first <- meta$Sequencing.run == "first"
second <- meta$Sequencing.run == "second"
lvs <- c(paste(meta$Object[first],
               meta$Sampling.number[first])[order(as.numeric(meta$Sampling.number[first]))],
         paste(meta$Object[second],
               meta$Sampling.number[second])[order(meta$Sampling.number[second])])
aframe$label <- factor(aframe$label, levels = lvs)

aframe$Institution[aframe$Collection == "RF"] <- "RF"
aframe$Institution <- factor(aframe$Institution, levels = c("MfN", "SPK", "RF"))

ggplot(aframe,
       aes(x = label, y = value, fill = variable)) + 
  facet_grid(~ interaction(Sequencing.run, Institution), 
             scales = "free_x", space = "free_x") +
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_manual(values = farben, name = "Phylum") +
  ylab("Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../../figs/phylum_abundance_barplot_all_samples_after_filtering.pdf",
       height = 8,
       width = 12)

# Discover touch signature in MfN data ####
tmp <- meta[which(meta$Institution == "MfN" &
                    meta$Sequencing.run == "first" &
                    meta$Comments != "reference" &
                    meta$total_reads > 500),]

ok_samples <- tmp$sample.id
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros >= 2)

normalized <- t(t(counts_tmp[rows_ok, ])/colSums(counts_tmp[rows_ok, ]))
normalized <- sqrt(normalized)

pca <- prcomp(t(normalized))
var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 3)
aframe <- data.frame(pca$x,
                     tmp[colnames(normalized), ])

aframe$group <- aframe$Collection
aframe$group[-grep("untouched", aframe$Comments)] <- "touched"

ggplot(aframe, aes(PC1, PC2,
                   color = group, group = group,
                   shape = Collection)) +
  labs(title = "MfN",
       subtitle = "Touched molluscs and fossils cluster together") +
  xlab(paste("PC1", var_explained[1], "%")) +
  ylab(paste("PC2", var_explained[2], "%")) +
  geom_point() +
  theme_classic()
ggsave("../../figs/molluscs_fossils_pca.pdf",
       height = 4,
       width = 5)

# Run DE - touched vs untouched molluscs & fossils ####
library(DESeq2)
aframe$touched <- !is.na(aframe$Touched)
des <- DESeqDataSetFromMatrix(counts_tmp[rowSums(counts_tmp) > 0, ] + 1,
                              colData = aframe,
                              design = ~ Collection + touched)
des <- DESeq(des)
res <- results(des)
res <- res[sort.list(res$pvalue), ]
aframe <- data.frame(res, TAX[rownames(res), ])

ggplot(aframe, aes(log2FoldChange, -log10(pvalue),
                   color = Genus == "g__Corynebacterium")) +
  geom_point() +
  theme_classic()

res_deseq <- res

# Run DE - touched vs untouched molluscs & fossils ####
res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x, meta[colnames(normalized), c("Comments", "Object", "Collection")])
  aframe$touched <-"yes"
  aframe$touched[grep("untouched", aframe$Comments)] <- "no"
  afit <- try(summary(lm(expr ~ Collection + touched, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  coefficients(afit)[3, c(1, 4)]
}))
colnames(res) <- c("coef", "p_value")
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res[, 2]),]
res$p_adjusted <- p.adjust(res$p_value, method = "BH")

sig <- rownames(res)[res$p_adjusted < 0.25 & res$coef > 0]
meta$touch_signature <- colSums(counts[sig,])/meta$total_reads

anno_col <- meta[colnames(normalized), c("Comments", "Collection", "Object", "total_reads", "touch_signature")]
anno_col$touched <-"yes"
anno_col$touched[grep("untouched", anno_col$Comments)] <- "no"
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(tax_info[rownames(normalized),])
rownames(anno_row) <- rownames(normalized)

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollucs` = "blue", `VertPal` = "darkblue")
object <- c(`Chiton` = "bisque3", `Tendaguru` = "darksalmon", `Triton horn`= "yellow", `Mobbi collection`= "brown")
ann_colors <- list(touched = touched,
                   collection = collection,
                   object = object)

p <- pheatmap(normalized[sig, order(anno_col$Collection, anno_col$touched, anno_col$Object)],
         cluster_cols = F,
         gaps_col = table(anno_col$Collection)[1],
         show_rownames = F, show_colnames = T,
         annotation_col = anno_col[, c("Comments", "Collection", "Object")],
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")],
         annotation_colors = ann_colors)
ggsave(p, filename = "../../figs/heatmap_molluscs_fossils_touched.pdf",
       width = 9,
       height = 12)

# Compare DESeq2 vs differential abundance ####
ok <- intersect(rownames(res_deseq),
                rownames(res))
aframe <- data.frame(res[ok,], res_deseq[ok,])
ggplot(aframe, aes(coef, log2FoldChange)) +
  labs(y = "Log2 fold change [DESeq2]",
       x = "Differential abundance coefficient") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  ggpubr::stat_cor() +
  theme_classic()
ggsave(filename = "../../figs/diff_abundance_vs_deseq2.pdf",
       width = 8,
       height = 8)

# MfN molluscs & fossils - Plot single ASV feature ####
aframe <- data.frame(anno_col,
                     cutibacterium = normalized["5b9edeb0a187d9a0ce1edf9c60bbb7b5",])

ggplot(aframe, aes(touched, cutibacterium, color = touched)) +
  facet_wrap(~ Collection) +
  labs(y = "Relative abundance Corynebacterium") +
  geom_boxplot() + geom_point() +
  theme_classic()
ggsave("../../figs/Corynebacterium_fossils_molluscs_touched.pdf",
       height = 5,
       width = 4)

# MfN molluscs & fossils - plot touch signature ####
nom <- c(Mollucs = "Molluscs", VertPal = "Fossils")
aframe$collection <- nom[aframe$Collection]
ggplot(aframe, aes(touched, touch_signature, color = touched)) +
  facet_wrap(~ collection) +
  labs(y = "Touch signature") +
  geom_boxplot() + geom_point() +
  theme_classic()
ggsave("../../figs/Touch_signature_fossils_molluscs_touched.pdf",
       height = 5,
       width = 4)

# Check if fingerprint is lost ####
ok_samples <- rownames(aframe)[aframe$touched == "no"]
counts_tmp <- counts[, ok_samples]
zeros <- apply(counts_tmp, 1, function(x) sum(x > 0))
rows_ok <- which(zeros >= 2)

normalized <- t(t(counts_tmp[rows_ok, ])/colSums(counts_tmp[rows_ok, ]))
normalized <- sqrt(normalized)

res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x, meta[colnames(normalized), c("Comments", "Object", "Collection")])
  afit <- try(summary(lm(expr ~ Collection, data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  coefficients(afit)[2, c(1, 4)]
}))
colnames(res) <- c("coef", "p_value")
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res[, 2]),]
res$p_adjusted <- p.adjust(res$p_value, method = "BH")

res$label <- res$Genus
res$label <- gsub("g__", "", fixed = T, res$label)

res$Phylum <- gsub("p__", "", fixed = T, res$Phylum)

res$label[11:nrow(res)] <- NA

ggplot(res, aes(coef, -log10(p_value), label = label, color = Phylum)) +
  geom_point() + ggrepel::geom_label_repel() +
  labs(y = "-log10 p-value",
       x = "Coefficient") +
  theme_classic()

plot_ASV <- function(asv){
  aframe <- data.frame(expr = normalized[asv, ],
                       meta[colnames(normalized), c("Comments", "Object", "Collection")])
  
  ggplot(aframe, aes(Collection, expr, color = Object)) +
    labs(title = res[asv, ]$label,
         y = "Relative abundance",
         x = "Touched by") +
    geom_boxplot() + geom_point() +
    theme_classic()  
}

# Define functions ####
get_subset <- function(samples){
  ok_samples_2nd <- samples
  
  counts_tmp <- counts[, ok_samples_2nd]
  counts_tmp <- counts_tmp[which(rowSums(counts_tmp) > 10), ]
  
  ok_samples_2nd <- colnames(counts_tmp)
  
  normalized <- t(t(counts_tmp)/colSums(counts_tmp))
  normalized <- sqrt(normalized)
  
  pca <- prcomp(t(normalized), scale. = T)
  var_explained <- signif(((pca$sdev/sum(pca$sdev))[1:2])*100, 2)
  
  aframe_tmp <- data.frame(pca$x[, 1:5],
                           meta[ok_samples_2nd, ])
  
  list(aframe_tmp, var_explained[1:2])
}

run_pca <- function(samples, color = "Object", label = "sample.id"){
  
  aframe_tmp <- get_subset(samples = samples)[[1]]
  var_explained <- get_subset(samples = samples)[[2]]
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
  counts <- counts[which(apply(counts, 1, function(x) sum(x > 0)) >= 2), ]
  sqrt(t(t(counts)/colSums(counts))) 
}

# SPK analysis - height profile ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Collection == "VAM" & 
                   meta$Sequencing.run == "second" &
                   meta$Material == "Glazed Ceramic" &
                   is.na(meta$Comments)), "sample.id"]

aframe_tmp <- get_subset(samples = ok)[[1]]
correl <- cor.test(aframe_tmp$PC1, as.numeric(aframe_tmp$Height))
p1 <- ggplot(aframe_tmp,
       aes(PC1, as.numeric(Height), color = Height, label = sample.id)) +
  geom_point() + geom_label() +
  labs(title = "SPK: exhibition hall",
       subtitle = paste0("Pearson Rho: ", signif(correl$estimate, 2),
                         "P: ", signif(correl$p.value, 2)),
       y = "Height [m]",
       x = "PC1") +
  ggpubr::stat_cor() +
  theme_classic()
summary(lm(PC1 ~ log10(total_reads) + as.numeric(Height), data = aframe_tmp))

correl <- cor.test(aframe_tmp$PC1, as.numeric(aframe_tmp$touch_signature))
p2 <- ggplot(aframe_tmp,
       aes(touch_signature, as.numeric(Height), color = Height, label = sample.id)) +
  geom_point() + geom_label() +
  labs(title = "SPK exhibition hall",
       subtitle = paste0("Pearson Rho: ", signif(correl$estimate, 2),
                         "P: ", signif(correl$p.value, 2)),
       y = "Height [m]",
       x = "Touch signature") +
  ggpubr::stat_cor() +
  theme_classic()
summary(lm(touch_signature ~ as.numeric(Height), data = aframe_tmp))

p <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(p, filename = "../../figs/Height_profile_vs_PC1_and_touch_signature.pdf",
       height = 5,
       width = 11)

# Supplemental figure ####
nom <- c(Stier = "Bull", Drache = "Dragon")
aframe_tmp$label <- nom[aframe_tmp$Object]
ggplot(aframe_tmp, aes(PC1, PC2, color = Color, label = label)) +
  geom_point() + ggrepel::geom_label_repel() +
  labs(title = "Pergamon museum: Ishtar gate",
       subtitle = "Glazed ceramic tiles") +
  scale_color_manual(values = c("#0D1C3D", "#956427")) +
  theme_classic()
ggsave(filename = "../../figs/pca_yellow_blue_tiles.pdf",
       height = 6,
       width = 7)

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
res_height <- data.frame(res_height,
                         tax_info[rownames(res_height), ])

anno_col <- meta[colnames(normalized), ]
anno_col$Height <- as.numeric(anno_col$Height)
rownames(anno_col) <- colnames(normalized)

anno_row <- data.frame(tax_info[rownames(normalized),])
rownames(anno_row) <- rownames(normalized)

sig <- rownames(res_height)[1:30]

p <- pheatmap(normalized[sig, order(anno_col$Height)],
         cluster_cols = F,
         scale = "row", breaks = seq(-1.5, 1.5, length = 100),
         show_rownames = F, show_colnames = F,
         annotation_col = anno_col[, c("Height", "Object", "Spot")],
         color = colorRampPalette(c("grey", "red"))(101),
         annotation_row = anno_row[, c("Phylum", "Genus")])
ggsave(p, filename = "../../figs/heatmap_height_profile.pdf",
       width = 9,
       height = 12)

# SPK analysis - blue vs yellow tiles ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Type.of.room == "exhibition" &
                   meta$Material == "Glazed Ceramic" &
                   is.na(meta$Comments)), "sample.id"]
aframe_tmp <- get_subset(samples = ok)[[1]]
run_pca(samples = ok, color = "Touched", label = "Object")
run_pca(samples = ok, color = "Color", label = "Object")

aframe_tmp <- get_subset(samples = ok)[[1]]
ggplot(aframe_tmp, aes(PC1, touch_signature)) +
  labs(y = "Touch signature",
       x = "PC1") +
  geom_point() +
  theme_classic()
summary(lm(PC1 ~ is.na(Touched) + Color, data = aframe_tmp))

# SPK analysis - Confirm touch signature on SPK gate lions ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   is.na(meta$Comments.2) &
                   meta$Material == "Stone"), "sample.id"]

aframe_tmp <- get_subset(samples = ok)[[1]]

aframe_tmp$group <- "small Gate lion" 
aframe_tmp$group[grep("Sam'al", aframe_tmp$Object)] <- "large Gate lion"

aframe_tmp$side <- "left" 
aframe_tmp$side[grep("right", aframe_tmp$Object)] <- "right"

afit <- summary(lm(touch_signature ~ group + is.na(Touched), data = aframe_tmp))
pval <- signif(coefficients(afit)[3, 4], 2)

ggplot(aframe_tmp,
       aes(Spot, touch_signature,
           color = Spot,
           group = Spot, shape = side)) +
  labs(title = "Gate lions",
       subtitle = paste("Multiple regression P", pval)) +
  facet_wrap(~ group, scales = "free_x") +
  geom_boxplot() +
  geom_point() +
  theme_classic()
ggsave("../../figs/gate_lions_touch_signature.pdf",
       height = 5,
       width = 5)

# SPK analysis - STS vs RRP ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Comments %in% c("StS", "visitors and StS", "RRP")), "sample.id"]

aframe_tmp <- get_subset(samples = ok)[[1]]
var_explained <- get_subset(samples = ok)[[2]]

p1 <- ggplot(aframe_tmp, aes(PC1, PC2,
                       shape = Material,
                       color = Comments,
                       label = Object)) +
  geom_point() + ggrepel::geom_label_repel() +
  xlab(paste("PC1", var_explained[1], "%")) + 
  ylab(paste("PC2", var_explained[2], "%")) + 
  theme_classic()

aframe_tmp$Material <- factor(aframe_tmp$Material,
                              levels = c("Polymer", "Glazed Ceramic", "Stone"))

p2 <- ggplot(aframe_tmp, aes(Material, PC2, color = Material)) +
  geom_boxplot() + geom_point() +
  ggpubr::stat_compare_means(comparisons = list(c(1, 2),
                                                c(1, 3)),
                             method = "t.test") +
  scale_color_manual(values = c("brown", "orange", "red")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p3 <- ggplot(aframe_tmp, aes(Comments, PC1, color = Comments)) +
  geom_boxplot() + geom_point() +
  ggpubr::stat_compare_means(comparisons = list(c(1, 2)),
                             method = "t.test") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

lay <- rbind(c(1,1,2,3))
p <- gridExtra::grid.arrange(p1, p3, p2, nrow = 1, layout_matrix = lay)

ggsave(p, filename = "../../figs/sts_vs_rrp.pdf",
       height = 6, width = 14)

# Run DE - RRP vs StS ####
counts_tmp <- counts[, ok]
counts_tmp <- counts_tmp[which(rowSums(counts_tmp) > 10), ]

ok_samples_2nd <- colnames(counts_tmp)

normalized <- t(t(counts_tmp)/colSums(counts_tmp))
normalized <- sqrt(normalized)

res <- t(apply(normalized, 1, function(x){
  aframe <- data.frame(expr = x,
                       meta[colnames(normalized), c("Comments", "Object", "Material")])
  
  afit <- try(summary(lm(expr ~ Material + (Comments == "RRP"), data = aframe)))
  if(class(afit) == "try-error") return(c(NA, NA))
  coefficients(afit)[4, c(1, 4)]
}))
colnames(res) <- c("coef", "p_value")
res <- data.frame(res, TAX@.Data[rownames(res),])
res <- res[sort.list(res[, 2]),]
res$p_adjusted <- p.adjust(res$p_value, method = "BH")

res$label <- res$Genus
res$label <- gsub("g__", "", fixed = T, res$label)

res$Phylum <- gsub("p__", "", fixed = T, res$Phylum)

res$label[11:nrow(res)] <- NA

p0 <- ggplot(res, aes(coef, -log10(p_value), label = label, color = Phylum)) +
  geom_point() + ggrepel::geom_label_repel() +
  labs(y = "-log10 p-value",
       x = "Coefficient") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73",
                                "#0072B2", "#CC79A7")) +
  theme_classic()

# Plot ASV with differential abundance bw subjects ####
plot_ASV <- function(asv){
  aframe <- data.frame(expr = normalized[asv, ],
                       meta[colnames(normalized), c("Comments", "Object", "Material")])
  
  aframe <- aframe[aframe$Comments %in% c("StS", "RRP"), ]
  aframe$Object <- gsub(" RRP", "", fixed = T, aframe$Object)
  aframe$Object <- gsub(" StS", "", fixed = T, aframe$Object)
  aframe$Object <- gsub(" right Sam'al", "", fixed = T, aframe$Object)
  aframe$Object <- gsub(" left Sam'al", "", fixed = T, aframe$Object)
  
  ggplot(aframe, aes(Comments, expr, color = Comments)) +
    facet_wrap(~ Material) +
    labs(title = res[asv, ]$label,
         y = "Relative abundance",
         x = "Touched by") +
    geom_boxplot() + geom_point() +
    scale_color_manual(values = c("#F3766E", "#2AB34B")) +
    theme_classic()  
}
p1 <- plot_ASV(asv = "08797478fa4139b217d0c5cc2d425501")
p2 <- plot_ASV(asv = "f13e5db13ee4930ce7fd9d44939141b9")

p <- grid.arrange(p0, p1, p2, nrow = 1, layout_matrix = lay)
ggsave(p, filename = "../../figs/example_differential_sts_vs_rrp.pdf",
       height = 8, width = 16)

# SPK analysis - WIP ####
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Material == "Glazed Ceramic"), "sample.id"]
run_pca(samples = ok, color = "Type.of.room", label = "Object")
run_pca(samples = ok, color = "Touched", label = "Object")
ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "second" &
                   meta$Material == "Glazed Ceramic"), "sample.id"]
ok <- ok[setdiff(1:length(ok),
                 which(meta[ok, "Comments"] %in% c("StS", "RRP")))]
run_pca(samples = ok, color = "Type.of.room", label = "Spot")
run_pca(samples = ok, color = "touch_signature", label = "Spot")

ok <- meta[which(meta$Institution == "SPK" &
                   meta$Sequencing.run == "first" &
                   meta$Object == "Zeus Sosipolis Tempel, Hellenistic Hall"), "sample.id"]

aframe_tmp <- get_subset(samples = ok)[[1]]
ggplot(aframe_tmp, aes(Spot, touch_signature, color = Touched)) +
  labs(title = "Zeus Sosipolis Tempel") +
  geom_boxplot() + geom_point() +
  theme_classic()
