def align_inputs(wildcards):

    fastq_path = path.join(fastq_dir, "{sample}_{{pair}}.fastq.gz".format(
        sample = wildcards.sample))
    pairs = ["R1", "R2"] if is_paired else [""]
    return expand(fastq_path, pair=pairs)


if is_no_index:
    ruleorder: bowtie_index > MapSplice

    rule bowtie_index:
        output:
            temp(touch(path.join(fusion_dir, '.build_index.done')))
        params:
            ref_fa = config['bowtie']['ref_fa'],
            label = config['bowtie']['name']
        threads:
            config['bowtie']['threads']
        conda:
            path.join(pipe_path, 'envs', 'bowtie.yaml')
        log:
            path.join(log_dir, 'bowtie_index', 'bowtie.log')
        shell:
            'bowtie-build -f {params.ref_fa} {params.label} '
            '--threads {threads} 2> {log}'


    rule MapSplice:
        input:
            fq = align_inputs,
            tmp = path.join(fusion_dir, '.build_index.done')
        output:
            fq_1 = temp(path.join(fastq_dir, '{sample}_R1.fastq')),
            fq_2 = temp(path.join(fastq_dir, '{sample}_R2.fastq')),
            fus = path.join(fusion_dir, '{sample}', 'fusions_well_annotated.txt')
        params:
            reflib = config['MapSplice']['reflib'],
            bw = config['bowtie']['name'],
            gtf = config['MapSplice']['gtf'],
            out_dir = path.join(fusion_dir, '{sample}'),
            options = format_options(config['MapSplice']['options'])
        threads:
            config['MapSplice']['threads']
        conda:
            path.join(pipe_path, 'envs', 'MapSplice.yaml')
        log:
            path.join(log_dir, 'MapSplice', '{sample}.log')
        shell:
            'zcat {input.fq[0]} > {output.fq_1} && '
            'zcat {input.fq[1]} > {output.fq_2} && '
            'mapsplice.py {params.options} -c {params.reflib} -x {params.bw} '
            '-o {params.out_dir} -p {threads} --gene-gtf {params.gtf} '
            '-1 {output.fq_1} -2 {outnput.fq_2} 2> {log}'


else:

    rule MapSplice:
        input:
            fq = align_inputs
        output:
            fq_1 = temp(path.join(fastq_dir, '{sample}_R1.fastq')),
            fq_2 = temp(path.join(fastq_dir, '{sample}_R2.fastq')),
            fus = path.join(fusion_dir, '{sample}', 'fusions_well_annotated.txt'),
            bam = path.join(fusion_dir, '{sample}', 'alignments.bam')
        params:
            reflib = config['MapSplice']['reflib'],
            bw = config['bowtie']['name'],
            gtf = config['MapSplice']['gtf'],
            out_dir = path.join(fusion_dir, '{sample}'),
            options = format_options(config['MapSplice']['options'])
        threads:
            config['MapSplice']['threads']
        conda:
            path.join(pipe_path, 'envs', 'MapSplice.yaml')
        log:
            path.join(log_dir, 'MapSplice', '{sample}.log')
        shell:
            'zcat {input.fq[0]} > {output.fq_1} && '
            'zcat {input.fq[1]} > {output.fq_2} && '
            'mapsplice.py {params.options} -c {params.reflib} -x {params.bw} '
            '-o {params.out_dir} -p {threads} --gene-gtf {params.gtf} '
            '-1 {output.fq_1} -2 {output.fq_2} 2> {log}'


rule bam_sort_index:
    input:
        path.join(fusion_dir, '{sample}', 'alignments.bam')
    output:
        path.join(fusion_dir, '{sample}', 'alignments.sort.bam')
    threads: 
        config['bam_sort_index']['threads']
    shell:
        'samtools sort -o {output} -@ {threads} {input} && '
        'samtools index {output}'

