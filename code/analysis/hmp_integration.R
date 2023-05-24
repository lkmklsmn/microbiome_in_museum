# Load R libs ####
library(qiime2R)
library(phyloseq)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)

# Load HMP data ####
otus <- read.table("../../data/otu_table_psn_v35.txt",
                   sep = '\t',
                   header = T,
                   row.names = 1,
                   check = F,
                   skip = 1,
                   comment = '')

tax <- otus[, ncol(otus)]
tax <- do.call(rbind, lapply(tax, function(x) strsplit(x, ";", fixed = T)[[1]][1:6]))
otus <- data.matrix(otus[, -ncol(otus)])

meta_hmp <- read.delim("../../data/v13_map_uniquebyPSN.txt")
rownames(meta_hmp) <- meta_hmp$X.SampleID
otus <- otus[, match(meta_hmp$X.SampleID, colnames(otus))]

meta_hmp$location <- "oral cavity"
meta_hmp$location[meta_hmp$HMPbodysubsite %in% c("Anterior_nares")] <- "nasal cavity"
meta_hmp$location[meta_hmp$HMPbodysubsite %in% c("Right_Retroauricular_crease",
                                             "Left_Retroauricular_crease",
                                             "Right_Antecubital_fossa",
                                             "Left_Antecubital_fossa")] <- "skin"
meta_hmp$location[meta_hmp$HMPbodysubsite %in% c("Stool")] <- "gastrointestinal tract"
meta_hmp$location[meta_hmp$HMPbodysubsite %in% c("Mid_vagina",
                                             "Posterior_fornix",
                                             "Vaginal_introitus")] <- "urogenital tract"

# Calculate proportions at genus level ####
asplit <- split(1:nrow(tax), tax[, 6])
asplit <- asplit[which(unlist(lapply(asplit, length)) > 1)]
norm <- do.call(rbind, lapply(asplit, function(x) colSums(otus[x, ])))
norm <- t(t(norm)/colSums(norm))
norm_hmp <- norm

# Load Lounsbery data ####
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

aframe <- data.frame(meta)
aframe$Institution[aframe$Collection == "RF"] <- "RF"
aframe$group <- paste(aframe$Institution)
aframe$group[aframe$group == "DSMZ"] <- "Controls"

aframe$group[aframe$Sampling.number %in% c("Pergamonmuseum  Negativprobe",
                                           "RF Negativprobe",
                                           "2VAM NC Werkstatt")] <- "Controls"

# Remove sample w low counts or ASV features detected ####
bad <- which(aframe$num_detected < 100 | aframe$total_reads < 500 | aframe$group == "Controls")
meta <- meta[-bad, ]
counts <- counts[, -bad]

# Normalize Lounsbery data ####
asplit <- split(rownames(tax_info), tax_info[, "Genus"])
asplit <- asplit[which(unlist(lapply(asplit, length)) > 1)]
norm <- do.call(rbind, lapply(asplit, function(x) colSums(counts[x, ])))
norm <- t(t(norm)/colSums(norm))
norm_lounsbery <- norm


# Plot single example ###
ok <- intersect(rownames(norm_hmp), rownames(norm_lounsbery))
asplit <- split(1:ncol(norm_hmp), meta_hmp$location)
norm_hmp_mean <- do.call(cbind, lapply(asplit, function(x) rowMeans(norm_hmp[, x], na.rm = T)))

aframe <- data.frame(touched = norm_lounsbery[ok, "L50"],
                     skin = norm_hmp_mean[ok, "skin"], label = ok)
aframe$label[aframe$touched < 0.06] <- NA
aframe$label <- gsub("g__", "", fixed = T, aframe$label)
ggplot(aframe, aes(touched, skin, label = label)) +
  labs(x = "Touched mollusc [MfN]",
       y = "Average skin profile [HMP]") +
  geom_point() + ggrepel::geom_label_repel() +
  ggpubr::stat_cor() +
  theme_classic()
ggsave(filename = "../../figs/scatter_plot_hmp_touched_mollusc.pdf",
       width = 5,
       height = 5)

# Heatmap avg levels ####
ok <- intersect(rownames(norm_hmp), rownames(norm_lounsbery))
correl <- cor(norm_hmp[ok,], norm_lounsbery[ok, ],
              use = 'pairwise.complete')

asplit <- split(1:ncol(correl), paste(meta$Collection, meta$Sequencing.run, is.na(meta$Touched)))
asplit <- asplit[c("Mollucs first TRUE", "Mollucs first FALSE",
                   "VertPal first TRUE", "VertPal first FALSE")]

correl_means <- do.call(cbind, lapply(asplit, function(x) rowMeans(correl[, x], na.rm = T)))

asplit <- split(1:nrow(correl), meta_hmp$location)
correl_means <- do.call(rbind, lapply(asplit, function(x) colMeans(correl_means[x, ], na.rm = T)))

anno_col <- data.frame(collection = c("Mollusc", "Mollusc", "Fossil", "Fossil"),
                       touched = c("no", "yes", "no", "yes"),
                       row.names = colnames(correl_means))

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollusc` = "blue", `Fossil` = "darkblue")
ann_colors <- list(touched = touched,
                   collection = collection)

p <- pheatmap(correl_means,
         annotation_col = anno_col,
         color = viridisLite::inferno(20),
         display_numbers = T,
         show_colnames = F,
         annotation_colors = ann_colors)

ggsave(p, filename = "../../figs/heatmap_avg_hmp.pdf",
       width = 9,
       height = 12)

# Heatmap sample level ####
ok <- intersect(rownames(norm_hmp), rownames(norm_lounsbery))
correl <- cor(norm_hmp[ok,], norm_lounsbery[ok, ],
              use = 'pairwise.complete')

asplit <- split(1:nrow(correl), meta_hmp$location)
correl_means <- do.call(rbind, lapply(asplit, function(x) colMeans(correl[x, ], na.rm = T)))

ok <- meta$sample.id[which(meta$Collection %in% c("VertPal", "Mollucs") &
                             meta$Sequencing.run == "first")]
aframe <- meta[ok,]
aframe$touched <- "no"
aframe$touched[!is.na(aframe$Touched)] <- "yes"
ok <- aframe$sample.id[order(aframe$Collection, aframe$touched)]

anno_col <- aframe
anno_col$Object[grep("Mobbi", anno_col$Object)] <- "Mobbi collection"

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollucs` = "blue", `VertPal` = "darkblue")
object <- c(`Chiton` = "bisque3", `Tendaguru` = "darksalmon", `Triton horn`= "yellow", `Mobbi collection`= "brown")
ann_colors <- list(touched = touched,
                   Collection = collection,
                   Object = object)

p <- pheatmap(correl_means[, ok],
              annotation_col = anno_col[, c("touched", "Object", "Collection")],
              #display_numbers = T,
              cluster_cols = F,
              gaps_col = 6,
              color = viridisLite::inferno(50),
              annotation_colors = ann_colors,
              show_colnames = F)

ggsave(p, filename = "../../figs/heatmap_hmp.pdf",
       width = 9,
       height = 6)

# Heatmap total ####
ok <- intersect(rownames(norm_hmp), rownames(norm_lounsbery))
correl <- cor(norm_hmp[ok,], norm_lounsbery[ok, ],
              use = 'pairwise.complete')

asplit <- split(1:nrow(meta_hmp), meta_hmp$location)
samples <- unlist(lapply(asplit, function(x) sample(x, 100)))

ok <- meta$sample.id[which(meta$Collection %in% c("VertPal", "Mollucs") &
                             meta$Sequencing.run == "first")]
aframe <- meta[ok,]
aframe$touched <- "no"
aframe$touched[!is.na(aframe$Touched)] <- "yes"
ok <- aframe$sample.id[order(aframe$Collection, aframe$touched)]

anno_col <- aframe
anno_col$Object[grep("Mobbi", anno_col$Object)] <- "Mobbi collection"

touched <- c(`no` = "#F8766D", `yes` = "#00BFC4")
collection <- c(`Mollucs` = "blue", `VertPal` = "darkblue")
object <- c(`Chiton` = "bisque3", `Tendaguru` = "darksalmon", `Triton horn`= "yellow", `Mobbi collection`= "brown")
ann_colors <- list(touched = touched,
                   Collection = collection,
                   Object = object)

p <- pheatmap(na.omit(correl[samples, ok]),
         cluster_cols = F,
         annotation_col = anno_col[, c("touched", "Object", "Collection")],
         annotation_colors = ann_colors,
         gaps_col = 6,
         color = viridisLite::inferno(50),
         annotation_row = meta_hmp[, c("HMPbodysubsite", "location")],
         show_rownames = F, show_colnames = F)

ggsave(p, filename = "../../figs/heatmap_hmp_total.pdf",
       width = 8,
       height = 10)
