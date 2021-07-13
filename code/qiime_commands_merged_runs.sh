cd /storage/thinc/data/projects/lounsbery/qiime_merged

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /storage/thinc/data/projects/lounsbery/qiime_merged/manifest.txt \
--output-path /storage/thinc/data/projects/lounsbery/qiime_merged/lounsbery.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime cutadapt trim-paired \
--i-demultiplexed-sequences /storage/thinc/data/projects/lounsbery/qiime_merged/lounsbery.qza \
--p-cores 8 \
# 341f          
--p-front-f CCTACGGGWGGCWGCAG \
# 803r
--p-front-r CTACCRGGGTATCTAATCC \
--o-trimmed-sequences primer-trimmed-VL_16S_PE.qza \
# Verbose to monitor trimming
--verbose \
# Save the output to a log
&> primer_trimming.log 

qiime demux summarize \
--i-data primer-trimmed-VL_16S_PE.qza \
--o-visualization primer-trimmed-VL_16S_PE.qzv \

qiime dada2 denoise-paired \
--p-n-threads 8 \
--i-demultiplexed-seqs primer-trimmed-VL_16S_PE.qza \
--p-trunc-len-f 220 \
--p-trunc-len-r 280 \
--output-dir DADA2_denoising_output \
--verbose \
&> DADA2_denoising.log
