# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(microbiome)
library(colorspace)

# Load new data ####
load("../../data/Touched_DADA2_objects_2022-07-28.RData")
counts_new <- data.matrix(seqtab.final[, rownames(metadata)])
rownames(counts_new) <- seqtab.final$Sequence

# Load old data ####
counts <- read_qza(
  "/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/table.qza")$data
counts_old <- counts

# Load metadata ####
meta <- readxl::read_excel("/Users/lukas/OneDrive/Miko/Lounsbery/data/both_runs_combined/sample_sheet_data_science.xlsx", sheet = "Tabelle1")
meta <- data.frame(meta)
rownames(meta) <- meta$sample.id

# Match count tables ####
colnames(counts_new) <- gsub("MfN_", "MfN", fixed = T, colnames(counts_new))
colnames(counts_new) <- gsub("LNC_", "LNC ", fixed = T, colnames(counts_new))
colnames(counts_new) <- gsub("L_", "L", fixed = T, colnames(counts_new))
colnames(counts_new) <- gsub("2MfN", "MfN", fixed = T, colnames(counts_new))
colnames(counts_new) <- gsub("_", " ", fixed = T, colnames(counts_new))

ok <- intersect(colnames(counts_new), colnames(counts_old))
counts_new <- counts_new[, ok]
counts_old <- counts_old[, ok]

# Plot ####
aframe <- data.frame(meta[match(ok, meta$sample.id), ],
                     total_old = colSums(counts_old),
                     total_new = colSums(counts_new),
                     detected_old = apply(counts_old, 2, function(x) sum(x > 0)),
                     detected_new = apply(counts_new, 2, function(x) sum(x > 0)))

p1 <- ggplot(aframe, aes(total_old, total_new,
                         color = Sequencing.run, label = sample.id)) +
  labs(title = "Total count",
       x = "Old",
       y = "New") +
  geom_point() + #ggrepel::geom_label_repel() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_classic()
                     
p2 <- ggplot(aframe, aes(detected_old, detected_new,
                         color = Sequencing.run, label = sample.id)) +
  labs(title = "ASVs detected",
       x = "Old",
       y = "New") +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 2)
