setwd("/storage/thinc/data/projects/lounsbery")
tmp <- list.files("Lounsbery@DSMZ/Seq_Data/seq2", full.names=T)
tmp <- tmp[grep("fastq.gz", tmp)]
dir <- getwd()

files <- unlist(lapply(tmp, function(x) paste(dir, x, sep = "/")))
sample <- gsub("Lounsbery@DSMZ/Seq_Data/seq2//", "", fixed = T, tmp)
sample <- sample[seq(1, length(sample), by = 2)]
sample <- gsub("_R1.fastq.gz", "", fixed = T, sample)

read_1 <- files[seq(1, length(files), by = 2)]
read_2 <- files[seq(2, length(files), by = 2)]

manifest <- data.frame(sample_id = sample,forward_absolute_filepath = read_1, reverse_absolute_filepath = read_2)
write.table(manifest, row.names = F, "qiime/manifest.txt", sep = "\t", quote = F)
