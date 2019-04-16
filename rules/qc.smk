if is_paired:

    rule fastqc:
        input:
            path.join(fastq_dir, '{sample}_{pair}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc', '{sample}_{pair}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc', '{sample}_{pair}_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


    rule multiqc:
        input:
            expand(path.join(qc_dir, 'fastqc', '{sample}_{pair}_fastqc.zip'),
                   sample=get_samples(), pair=['R1', 'R2'])
        output:
            path.join(qc_dir, 'multiqc_report.html')
        params:
            qc_dir=qc_dir,
            output_dir=qc_dir
        shell: 'multiqc --outdir {params.output_dir} {params.qc_dir}'

else:

    rule fastqc:
        input:
            path.join(fastq_dir, '{sample}.fastq.gz')
        output:
            path.join(qc_dir, 'fastqc', '{sample}_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


    rule multiqc:
        input:
            expand(path.join(qc_dir, 'fastqc', '{sample}_fastqc.zip'),
                   sample=get_samples())
        output:
            path.join(qc_dir, 'multiqc_report.html')
        params:
            qc_dir=qc_dir,
            output_dir=qc_dir
        shell: 'multiqc --outdir {params.output_dir} {params.qc_dir}'

