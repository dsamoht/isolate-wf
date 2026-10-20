# dsamoht/isolate-wf: Usage

## Introduction

isolate-wf builds a consensus assembly of a bacterial isolate from long reads with
[Autocycler](https://github.com/rrwick/Autocycler). Each sample's reads are QC'd, split into
subsets, and every subset is assembled with nine assemblers (Flye, Plassembler, Canu,
metaMDBG, miniasm, Myloasm, Raven, NECAT and NextDenovo). Autocycler then compresses,
clusters, trims and resolves those assemblies into a single consensus, which is reoriented
with Dnaapler and assessed with CheckM2.

## Samplesheet input

Pass a comma-separated file with four columns to `--input`.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

```csv title="samplesheet.csv"
sample_id,long_reads,short_reads_1,short_reads_2
test1,tests/data/test1_reads_ont_r10.fastq.gz,,
test2,tests/data/test2_reads_ont_r10.fastq.gz,,
```

| Column          | Description                                                                                           |
| --------------- | ----------------------------------------------------------------------------------------------------- |
| `sample_id`     | Sample identifier, unique across the sample sheet. Spaces are not allowed.                            |
| `long_reads`    | FastQ file of long reads, gzipped or not. Required.                                                   |
| `short_reads_1` | Optional FastQ file (gzipped or not) of forward paired-end reads, already QC'd.                       |
| `short_reads_2` | Optional FastQ file (gzipped or not) of reverse paired-end reads. Required if `short_reads_1` is set. |

The sample sheet is validated against [assets/schema_input.json](../assets/schema_input.json).
Rules that span several rows are checked before any task is submitted:

- Sample ids must be unique: they name every output directory of a sample.
- A read file cannot be shared by two samples.

> [!NOTE]
> Short read polishing is not yet implemented: the short read columns are validated but not used.

## Running the pipeline

```bash
nextflow run dsamoht/isolate-wf -profile docker --input ./samplesheet.csv --outdir ./results
```

Results, work directory and Nextflow cache are written to the current directory.

### Updating the pipeline

```bash
nextflow pull dsamoht/isolate-wf
```

### Reproducibility

Specify the pipeline version with `-r`, and record it in your methods:

```bash
nextflow run dsamoht/isolate-wf -r 0.1.0 -profile docker --input ./samplesheet.csv --outdir ./results
```

## Core Nextflow arguments

### `-profile`

Use one of the container profiles rather than installing the tools yourself. Profiles are
comma-separated and order matters, the last one wins.

- `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `wave`
- `conda` / `mamba` — only if none of the container engines are available, and only together with
  `--skip_checkm2`: CheckM2 needs a DIAMOND reference database that the pipeline supplies through
  its container image, and the Bioconda package ships none. Plassembler likewise needs its database,
  which only the `dsamoht/autocycler` image ships, so the Plassembler assemblies fail under conda
  (Autocycler carries on with the other assemblers)
- `drac` — slurm executor and resources for Digital Research Alliance of Canada clusters
- `test` — a minimal dataset that exercises the whole pipeline

### Resource limits

Memory and time requests grow with each retry, and are capped so that no task asks for more
than a machine has: 256 GB, 32 CPUs and 72 h by default, 249 GB, 64 CPUs and 168 h under `drac`.
Set `process.resourceLimits` in a config passed with `-c` to match your largest node.

### `-resume`

Restart a previous run, reusing the tasks whose inputs have not changed.

## Pipeline specific options

| Option                | Default   | Description                                                                             |
| --------------------- | --------- | --------------------------------------------------------------------------------------- |
| `--read_type`         | `ont_r10` | Read type passed to the assemblers: `ont_r9`, `ont_r10`, `pacbio_clr` or `pacbio_hifi`. |
| `--chopper_minlength` | `1000`    | Minimum read length kept by Chopper.                                                    |
| `--chopper_minq`      | `10`      | Minimum average read quality kept by Chopper.                                           |

Each assembler can be turned off: `--skip_flye`, `--skip_plassembler`, `--skip_canu`,
`--skip_metamdbg`, `--skip_miniasm`, `--skip_myloasm`, `--skip_raven`, `--skip_necat`,
`--skip_nextdenovo`. Every read subset is then assembled with the remaining ones only, and
at least one must stay on: skipping all of them stops the run at start-up. Fewer assemblers
make the run faster but give Autocycler less to build its consensus from.

Any step can be turned off: `--skip_qc`, `--skip_nanoplot`, `--skip_porechop`,
`--skip_dnaapler`, `--skip_checkm2`, `--skip_multiqc`. All default to `false`.

The full, always up-to-date list of parameters is generated from
[nextflow_schema.json](../nextflow_schema.json):

```bash
nextflow run dsamoht/isolate-wf --help                # summary
nextflow run dsamoht/isolate-wf --help_full           # every parameter, with help text
nextflow run dsamoht/isolate-wf --help --show_hidden  # include hidden parameters
nextflow run dsamoht/isolate-wf --help read_type      # help for a single parameter
nextflow run dsamoht/isolate-wf --version             # pipeline version
```
