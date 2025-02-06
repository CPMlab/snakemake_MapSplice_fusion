# RNA-seq Fusion Detection Pipeline (MapSplice)

## Overview
This pipeline automates RNA sequencing (RNA-seq) analysis for detecting fusion genes using MapSplice. It processes input FASTQ files, aligns reads, and identifies fusion transcripts in a structured workflow using Snakemake.

## Dependencies
The following software and Python libraries are required:

- **MapSplice**: Fusion detection tool for RNA-seq
- **Snakemake**: Workflow management system
- **Python**: Required for Snakemake execution
- **Pandas**: For sample data processing
- **NumPy**: For numerical operations
- **OS module**: For path handling

## Pipeline Structure
The pipeline consists of the following Snakemake rules:

### 1. `rule all`
Defines the final output, ensuring that fusion detection and alignments are successfully generated in `fusion_dir`.
- **Output Files**:
  - `fusions_well_annotated.txt`: Annotated fusion transcript report
  - `alignments.sort.bam`: Sorted alignment file

### 2. `get_samples()`
Extracts the list of sample names from `samples.tsv`, ensuring that all required samples are processed.

### 3. `format_options(options)`
Formats optional command-line arguments into a single string for MapSplice execution.

## Configuration
The `config.yaml` file should specify input paths, output directories, and analysis options. Example:

```yaml
samples: "samples.tsv"
options:
  paired: true
  no_index: false
path:
  data: "./data"
  pipeline: "./pipeline"
```

## Included Rules
The pipeline includes additional Snakemake rule files:
- `rules/qc.smk`: Quality control analysis
- `rules/MapSplice.smk`: RNA-seq alignment and fusion detection

## Execution
To run the pipeline, execute:
```bash
snakemake --cores <num_cores>
```
Replace `<num_cores>` with the number of available CPU cores.

## Output Files
- **Fusion Results**: `fusions_well_annotated.txt`
- **Aligned Reads**: `alignments.sort.bam`
- **Quality Control Report** (if enabled): `multiqc_report.html`
- **Log Files**: Stored in `log_dir`

## Troubleshooting
- Ensure `config.yaml` contains correct paths and options.
- Verify `samples.tsv` has the correct format with a `sample` column.
- Check `log_dir` for error messages.
- Ensure MapSplice is installed and properly configured.

## Author
This pipeline was developed for RNA-seq fusion detection using MapSplice and Snakemake.

For inquiries, please contact the developer or refer to the MapSplice documentation.

