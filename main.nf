nextflow.enable.dsl=2

/*
 End-to-end QC ∩ Demuxlet workflow for all pools.
 CSV input (--pools_csv) columns:
   pool,bam,barcodes,donor_vcf,qc_rds,outdir
 Outputs go to the per-pool outdir specified in pools.csv.
*/

params.pools_csv         = params.pools_csv         ?: 'pools.csv'
params.pop_vcf           = params.pop_vcf           ?: '/home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.gz'
params.ref_fasta         = params.ref_fasta         ?: '/home/pr422/RDS/live/Users/Parisa/refrence_genome/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa'
params.doublet_prior     = params.doublet_prior     ?: 0.05
params.geno_error_offset = params.geno_error_offset ?: 0.15
params.geno_error_coeff  = params.geno_error_coeff  ?: 0.0
params.alphas            = params.alphas            ?: '0,0.5'
params.pileup_script     = params.pileup_script     ?: "${projectDir}/scripts/01_run_pileup_1000G.sh"
params.demux_script      = params.demux_script      ?: "${projectDir}/scripts/02_run_demuxlet_sweep.sh"
params.qc_intersect_R    = params.qc_intersect_R    ?: "${projectDir}/scripts/qc_intersect.R"

workflow {
    // Resolve pools CSV and fail fast if missing
    def pools_file = file(params.pools_csv)
    if( !pools_file.exists() ) error "Pools CSV not found: ${pools_file}"

    // Samples channel as tuples (path objects for files)
    samples = Channel.fromPath(pools_file)
        .splitCsv(header:true)
        .map { row -> tuple(row.pool, file(row.bam), file(row.barcodes), file(row.donor_vcf), file(row.qc_rds), row.outdir) }

    // Broadcast pop_vcf/ref to every sample (skip pop normalization; use provided files)
    def pop_vcf = file(params.pop_vcf)
    def pop_vcf_idx = file(params.pop_vcf + ".tbi")
    def ref_fa  = file(params.ref_fasta)
    def ref_fai = file(params.ref_fasta + '.fai')
    sample_with_pop = samples.map { s -> tuple(s[0], s[1], s[2], s[3], s[4], s[5], pop_vcf, pop_vcf_idx, ref_fa, ref_fai) }

    donor_norm    = sample_with_pop | NORMALIZE_DONOR
    intersected   = donor_norm      | INTERSECT_SITES
    reheadered    = intersected     | REHEADER_SITES
    harmonized    = reheadered      | HARMONIZE_DONOR
    pileups       = harmonized      | PILEUP
    demux_calls   = pileups         | DEMUXLET
    qc_labels     = demux_calls     | QC_INTERSECT
    reports       = qc_labels       | REPORT

    reports.view()
}

process NORMALIZE_DONOR {
    tag { "donor_${pool}" }
    input:
      tuple val(pool), val(bam), val(barcodes), val(donor_vcf), val(qc_rds), val(outdir), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai)
    output:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path('donor.norm.bi.snp.vcf.gz')
    script:
    """
    bcftools +fixref ${donor_vcf} -- -f ${ref_fasta} -m flip \\
      | bcftools norm -f ${ref_fasta} -m-any \\
      | bcftools view -m2 -M2 -v snps \\
      | bgzip -c > donor.norm.bi.snp.vcf.gz
    tabix -f -p vcf donor.norm.bi.snp.vcf.gz
    """
}

process INTERSECT_SITES {
    tag { "intersect_${pool}" }
    input:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(donor_norm)
    output:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path('sites.intersect.vcf.gz'), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(donor_norm)
    script:
    """
    # Ensure donor VCF is indexed
    bcftools index -f ${donor_norm}
    bcftools isec -n=2 -w1 ${pop_vcf} ${donor_norm} -Oz -o sites.intersect.vcf.gz
    bcftools index -f sites.intersect.vcf.gz
    """
}

process REHEADER_SITES {
    tag { "reheader_${pool}" }
    input:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(sites_vcf), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(donor_norm)
    output:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path('sites.intersect.bamorder.vcf.gz'), path('sites.intersect.bamorder.vcf.gz.tbi'), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(donor_norm)
    script:
    """
    if [ ! -f "${ref_fai}" ]; then
      echo "[ERROR] Missing FAI index for reference: ${ref_fai}" >&2
      exit 1
    fi
    # Build fresh contig header lines in the exact FASTA (and BAM) order
    awk 'BEGIN{OFS=""} {printf("##contig=<ID=%s,length=%s>\\n",\$1,\$2)}' ${ref_fai} > contigs.txt
    # Strip existing contig lines from header
    bcftools view -h ${sites_vcf} | grep -v '^##contig=' > header.no_contig.txt
    # Insert contigs before the #CHROM line
    awk 'NR==FNR{c[++n]=\$0;next} /^#CHROM/{for(i=1;i<=n;i++) print c[i]} {print}' contigs.txt header.no_contig.txt > header.with_contig.txt
    # Reheader and convert back to bgzipped VCF
    bcftools reheader -h header.with_contig.txt ${sites_vcf} -o sites.intersect.bamorder.bcf
    bcftools view -Oz -o sites.intersect.bamorder.vcf.gz sites.intersect.bamorder.bcf
    bcftools index -t -f sites.intersect.bamorder.vcf.gz
    """
}

process HARMONIZE_DONOR {
    tag { "harmonize_${pool}" }
    input:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(sites_bamorder), path(sites_bamorder_idx), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(donor_norm)
    output:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(sites_bamorder), path(sites_bamorder_idx), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path('donor.harmonized.to_sites.vcf.gz')
    script:
    """
    # Ensure donor VCF is indexed
    bcftools index -f ${donor_norm}
    # Ensure sites VCF is indexed
    bcftools index -f ${sites_bamorder}
    bcftools isec -c all -n=2 -w1 ${donor_norm} ${sites_bamorder} -Oz -o donor.harmonized.to_sites.vcf.gz
    bcftools index -f donor.harmonized.to_sites.vcf.gz
    """
}

process PILEUP {
    tag { "pileup_${pool}" }
    input:
      tuple val(pool), val(bam), val(barcodes), val(qc_rds), val(outdir), path(sites_bamorder), path(sites_bamorder_idx), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(harmonized_donor)
    output:
      tuple val(pool), val(barcodes), val(qc_rds), val(outdir), path(sites_bamorder), path(sites_bamorder_idx), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(harmonized_donor), path('pileup_prefix.txt')
    script:
    """
    export OMP_NUM_THREADS=60
    prefix="${outdir}/${pool}_pileup_intersect"
    /bin/bash ${params.pileup_script} \
      --sam ${bam} \
      --barcodes ${barcodes} \
      --vcf ${sites_bamorder} \
      --out \$prefix
    echo \$prefix > pileup_prefix.txt
    """
}

process DEMUXLET {
    tag { "demux_${pool}" }
    input:
      tuple val(pool), val(barcodes), val(qc_rds), val(outdir), path(sites_bamorder), path(sites_bamorder_idx), path(pop_vcf), path(pop_vcf_index), path(ref_fasta), path(ref_fai), path(harmonized_donor), path(pileup_prefix_file)
    output:
      tuple val(pool), val(qc_rds), val(outdir), path('demuxlet_err015.best')
    script:
    """
    export OMP_NUM_THREADS=60
    PLP=\$(cat ${pileup_prefix_file})
    OUTDIR="${outdir}/demuxlet"
    mkdir -p "\$OUTDIR"
    /bin/bash ${params.demux_script} \
      --plp \$PLP \
      --vcf ${harmonized_donor} \
      --barcodes ${barcodes} \
      --outdir "\$OUTDIR" \
      --errs ${params.geno_error_offset} \
      --doublet-prior ${params.doublet_prior} \
      --alpha-list ${params.alphas}
    ln -sf "\$OUTDIR"/demuxlet_err${params.geno_error_offset}.best demuxlet_err015.best
    """
}

process QC_INTERSECT {
    tag { "qc_intersect_${pool}" }
    input:
      tuple val(pool), val(qc_rds), val(outdir), path(demux_best)
    output:
      tuple val(pool), val(qc_rds), val(outdir), path(demux_best), path("${pool}_final_labels.tsv")
    script:
    """
    Rscript ${params.qc_intersect_R} ${qc_rds} ${demux_best} ${pool}_final_labels.tsv
    mkdir -p ${outdir}
    cp ${pool}_final_labels.tsv ${outdir}/
    """
}

process REPORT {
    tag { "report_${pool}" }
    input:
      tuple val(pool), val(qc_rds), val(outdir), path(demux_best), path(final_labels)
    output:
      tuple path("${pool}_summary.tsv"), path("${pool}_summary.png")
    script:
    """
    Rscript ${projectDir}/scripts/report_counts.R \\
      ${pool} ${qc_rds} ${demux_best} ${final_labels} ${pool}_summary.tsv ${pool}_summary.png
    mkdir -p ${outdir}
    cp ${pool}_summary.tsv ${pool}_summary.png ${outdir}/
    """
}
