#!/usr/bin/bash
#SBATCH --mem=12000
#SBATCH --job-name qiime2
#SBATCH --ntasks=4
#SBATCH --output /storage/thinc/home/233454/slurm_out/qiime2-log-%J.out
#SBATCH --error /storage/thinc/home/233454/slurm_out/qiime2-log-%J.err
## get tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)

cd /storage/thinc/data/projects/lounsbery/both_runs_combined 

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest_file \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

qiime vsearch join-pairs \
--i-demultiplexed-seqs paired-end-demux.qza \
--o-joined-sequences joined-demux.qza

qiime demux summarize \
  --i-data joined-demux.qza \
  --o-visualization joined-demux.qzv


### 2)	QUALITY FILTERING

qiime quality-filter q-score \
--i-demux joined-demux.qza \
--o-filtered-sequences demux-joined-filtered.qza \
--o-filter-stats demux-joined-filter-stats.qza

qiime metadata tabulate \
  --m-input-file demux-joined-filter-stats.qza \
  --o-visualization demux-joined-filter-stats.qzv

qiime deblur denoise-16S \
--i-demultiplexed-seqs demux-joined-filtered.qza \
--p-min-size 1 \
--p-min-reads 2 \
--p-trim-length 160 \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza \
--p-sample-stats \
--p-jobs-to-start 2 \
--o-stats deblur-stats.qza

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv


### 3) CONSTRUCT PHYLOGENETIC TREE

qiime alignment mafft \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--p-n-threads 1

qiime feature-table tabulate-seqs \
  --i-data aligned-rep-seqs.qza \
  --o-visualization aligned-rep-seqs.qzv

qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

qiime feature-table tabulate-seqs \
--i-data masked-aligned-rep-seqs.qza \
--o-visualization masked-aligned-rep-seqs.qzv

qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--p-n-threads 1

qiime phylogeny midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza 


### 4)	TAXONOMIC CLASSIFICATION

qiime feature-classifier classify-sklearn \
--i-classifier ../qiime/classifiers/silva_132_99_v3v4_q2_2019-7.qza \
--p-n-jobs 1 \
--p-reads-per-batch 1000 \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza 

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv


