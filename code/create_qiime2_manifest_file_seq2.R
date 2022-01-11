setwd("/storage/thinc/data/projects/lounsbery")
tmp <- list.files("Lounsbery@DSMZ/Seq_Data_2/Sampling2_SeqData", full.names=T)
tmp <- tmp[grep("fastq.gz", tmp)]
dir <- getwd()

files <- unlist(lapply(tmp, function(x) paste(dir, x, sep = "/")))
sample <- gsub("Lounsbery@DSMZ/Seq_Data_2/Sampling2_SeqData/", "", fixed = T, tmp)
sample <- sample[seq(1, length(sample), by = 2)]
sample <- gsub("_R1.fastq.gz", "", fixed = T, sample)

read_1 <- files[seq(1, length(files), by = 2)]
read_2 <- files[seq(2, length(files), by = 2)]

manifest <- data.frame(sample_id = sample,forward_absolute_filepath = read_1, reverse_absolute_filepath = read_2)
write.table(manifest, row.names = F, "qiime/manifest_seq2.txt", sep = "\t", quote = F)
                       
# Fix all column names by switching, for example, sample_id to sample-id using vim or other editor
