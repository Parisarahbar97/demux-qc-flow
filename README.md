Demux + QC Intersection Pipeline 
=======================================

This project runs an end-to-end demultiplexing workflow for all pools, producing high-confidence singlet labels as the intersection of:
- popscle demuxlet (harmonized 1000G pileup + donor VCF, alpha 0/0.5, prior 0.05, ERR 0.15)
- scQC-flow transcriptomic QC (filters: `doublet_class == "singlet"`, `nuclear_fraction >= 0.4`, `percent.mt <= 20`)

Inputs (pools.csv)
------------------
`pool,bam,barcodes,donor_vcf,qc_rds,outdir` for each pool (S1A–S19B). Example paths used here:
- BAM: `/home/pr422/RDS/live/Users/Parisa/alex_output/epilep_cellranger_outputs/${POOL}_mapped/outs/possorted_genome_bam.bam`
- Barcodes (CellBender whitelist): `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest/${POOL}/${POOL}_cellbender_output/cellbender_out_cell_barcodes.csv`
- Donor VCF: `/home/pr422/RDS/live/Users/Parisa/alex_output/genotype_inputs2/${POOL}.vcf.gz`
- QC Seurat RDS: `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest/${POOL}/${POOL}_seurat_object.rds`
- Outdir: `/home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/demuxlet_results/${POOL}`

Global defaults (main.nf)
-------------------------
- Pop VCF: `/home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.gz`
- FASTA: `/home/pr422/RDS/live/Users/Parisa/refrence_genome/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa`
- Demuxlet: `alpha 0,0.5`, `--doublet-prior 0.05`, `--geno-error-offset 0.15`, `--geno-error-coeff 0.0`

Pipeline steps (per pool)
-------------------------
1. **Normalize donor VCF** – `bcftools +fixref|norm|view` vs `genome.fa` → `donor.norm.bi.snp.vcf.gz` (bgzip + tabix).  
2. **Intersect sites** – index donor VCF, `bcftools isec -n=2 -w1` against the population panel to get `sites.intersect.vcf.gz`.  
3. **Reheader to BAM/FASTA order** – rebuild contig header from `genome.fa.fai`, reheader/interconvert (BCF→VCF.gz) and index → `sites.intersect.bamorder.vcf.gz`.  
4. **Harmonize donor variants** – index donor/sites VCFs, `bcftools isec -c all -n=2 -w1` → `donor.harmonized.to_sites.vcf.gz` (indexed).  
5. **Pileup** – `scripts/01_run_pileup_1000G.sh` (popscle dsc-pileup, `OMP_NUM_THREADS=60`, CB/UB tags) produces pileup matrices and `pileup_prefix.txt`.  
6. **Demuxlet** – `scripts/02_run_demuxlet_sweep.sh` using the harmonized donor VCF + pileup prefix (`alpha 0/0.5`, `doublet-prior 0.05`, `geno-error-offset 0.15`) and writes `demuxlet_err015.best`.  
7. **QC intersection** – `scripts/qc_intersect.R` merges Seurat QC labels with demux calls to create `<pool>_final_labels.tsv` and copies it to the pool outdir.  
8. **Reporting** – `scripts/report_counts.R` summarises singlet/doublet counts into `<pool>_summary.tsv/.png`, again copied to the outdir.

Containers
----------
- QC_INTERSECT runs in `ghcr.io/johnsonlab-ic/sc_analysis:latest` (set in `nextflow.config`, docker enabled for local profile).  
- Pileup/demuxlet call their own Docker images inside the bash scripts (`parisa/demux:2.1`).

How to run (current production command)
---------------------------------------
The runs that are executing successfully use a scratch temp/work directory on `/home/pr422/RDS/ephemeral/parisa_tmp` and rely on Docker (local profile):

```bash
cd /home/pr422/RDS/live/Users/Parisa/demux-qc-flow
NXF_TEMP=/home/pr422/RDS/ephemeral/parisa_tmp/tmp \
nextflow run main.nf -profile local \
  --pools_csv /home/pr422/RDS/live/Users/Parisa/demux-qc-flow/pools.csv \
  -work-dir /home/pr422/RDS/ephemeral/parisa_tmp/work
```

Key notes:
- Always include `-resume` when restarting a run; all intermediate steps are cache-aware.
- The `pools.csv` shipped here has S3B commented/removed (no donor VCF available). Adjust as needed but keep the column order.
- Ensure `genome.fa.fai` is present alongside the FASTA (required by REHEADER_SITES).
- Outputs for each pool are written to the `outdir` path defined in `pools.csv` (e.g. `/home/pr422/RDS/live/Users/Parisa/EPILEP/demuxlet_results/S1A/`), containing pileup, demuxlet, `<pool>_final_labels.tsv`, and `<pool>_summary.{tsv,png}`.
