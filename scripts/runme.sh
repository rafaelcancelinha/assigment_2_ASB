#!/usr/bin/env bash
mkdir fastq_L1
mkdir fastq_L2_L3

python3 scripts/test_ncbi.py importante/ancona_list.txt 
mv *.fastq.gz fastq_L2_L3

python3 scripts/test_ncbi.py importante/napoles.txt
mv *.fastq.gz fastq_L2_L3

python3 scripts/test_ncbi.py importante/denamark.txt
mv *.fastq.gz fastq_L1


echo -e "sampleID\tforwardReads\treverseReads" > samplesheet_l1.tsv
for f1 in fastq_L1/*_1.fastq.gz; do
    sample=$(basename "$f1" _1.fastq.gz)
    f2="fastq_L1/${sample}_2.fastq.gz"
    if [ -f "$f1" ] && [ -f "$f2" ]; then
        echo -e "${sample}\t$(pwd)/${f1}\t$(pwd)/${f2}" >> samplesheet_l1.tsv
    fi
done



echo -e "sampleID\tforwardReads\treverseReads" > samplesheet_l23.tsv
for f1 in fastq_L2_L3/*_1.fastq.gz; do
    sample=$(basename "$f1" _1.fastq.gz)
    f2="fastq_L2_L3/${sample}_2.fastq.gz"
    if [ -f "$f1" ] && [ -f "$f2" ]; then
        echo -e "${sample}\t$(pwd)/${f1}\t$(pwd)/${f2}" >> samplesheet_l23.tsv
    fi
done


cat << 'INNER_EOF' > nextflow.config
params {
    trunclenf = 260
    trunclenr = 180
    max_ee = 3
    skip_dada_addspecies = true
    mergepairs_strategy = "merge"
    skip_cutadapt = true
}

process {
    cpus = 2
    memory = '6 GB'

    withName: 'NFCORE_AMPLISEQ:AMPLISEQ:DADA2_DENOISING' {
        cpus = 10
        memory = '12 GB'
    }

    withName: '.*DADA2_TAXONOMY.*' {
        cpus = 10
        memory = '12 GB'
    }
}
INNER_EOF


echo "=============  ARRANCAR PIPELINE L1 ============="
nextflow run nf-core/ampliseq \
  -profile conda \
  -c nextflow.config \
  --input ./samplesheet_l1.tsv \
  --metadata ./metadata_l1.tsv \
  --outdir ./resultados_L1 \
  --FW_primer GTGYCAGCMGCCGCGGTAA \
  --RV_primer GGACTACNVGGGTWTCTAAT
  -resume

echo "=============  ARRANCAR PIPELINE L2_L3 ============="
nextflow run nf-core/ampliseq \
  -profile conda \
  -c nextflow.config \
  --input ./samplesheet_l23.tsv \
  --metadata ./metadata_l23.tsv \
  --outdir ./resultados_L2_L3 \
  --FW_primer GTGYCAGCMGCCGCGGTAA \
  --RV_primer GGACTACNVGGGTWTCTAAT \
  -resume

echo "=============  TODOS OS PIPELINES TERMINADOS! ============="


Rscript scripts/cleaning_data.R 
Rscript scripts/analise_completa.R resultados_cleaning resultados_analise

