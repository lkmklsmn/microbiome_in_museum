#conda activate qiime2-2021.4
cd /storage/thinc/data/projects/lounsbery/qiime

qiime cutadapt trim-paired \
--i-demultiplexed-sequences lounsbery_seq2.qza \
--p-cores 16 \
--p-front-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--p-front-r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
--o-trimmed-sequences lounsbery_seq2_trimmed.qza \
--verbose \
&> primer_trimming.log 

qiime demux summarize \
--i-data lounsbery_seq2_trimmed.qza \
--o-visualization lounsbery_seq2_trimmed.qzv

qiime dada2 denoise-paired \
--p-n-threads 16 \
--i-demultiplexed-seqs lounsbery_seq2_trimmed.qza \
--p-trunc-len-f 285 \
--p-trunc-len-r 200 \
--output-dir DADA2_denoising_output \
--verbose \
&> DADA2_denoising.log

qiime metadata tabulate \
--m-input-file DADA2_denoising_output/denoising_stats.qza \
--o-visualization DADA2_denoising_output/denoising_stats.qzv

qiime feature-table tabulate-seqs \
--i-data DADA2_denoising_output/representative_sequences.qza \
--o-visualization DADA2_denoising_output/rep_seqs.qzv

qiime feature-table summarize \
--i-table DADA2_denoising_output/table.qza \
--o-visualization DADA2_denoising_output/table.qzv

#qiime feature-classifier classify-sklearn \
#--i-classifier ../genomes/rlappan/SILVA/classifier_16S_V3-V4.qza \
#--i-reads DADA2_denoising_output/representative_sequences.qza \
#--o-classification classified_rep_seqs.qza

# Download classifier
# wget -O "gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2018.2/common/gg-13-8-99-515-806-nb-classifier.qza"
# https://github.com/Jiung-Wen/q2-silva-V3V4classifier
qiime feature-classifier classify-sklearn \
--i-classifier classifiers/gg-13-8-99-515-806-nb-classifier.qza \
--i-reads DADA2_denoising_output/representative_sequences.qza \
--o-classification classified_rep_seqs.qza

qiime metadata tabulate \
--m-input-file classified_rep_seqs.qza \
--o-visualization classified_rep_seqs.qzv

