[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

```
 _           _       _                      __
(_)___  ___ | | __ _| |_ ___     __      __/ _|
| / __|/ _ \| |/ _` | __/ _ \____\ \ /\ / / |_
| \__ \ (_) | | (_| | ||  __/_____\ V  V /|  _|
|_|___/\___/|_|\__,_|\__\___|      \_/\_/ |_|

Nextflow wrapper for the Autocycler long-read consensus
genome assembly pipeline

     Github: https://github.com/dsamoht/isolate-wf
     Version: v0.1.0
```

## tl;dr

```bash
nextflow run main.nf --help # if working in the repo
```

or

```bash
nextflow run dsamoht/isolate-wf --help # nextflow pulls the repo automatically
```

## software dependencies

- Nextflow
- a container engine (i.e. docker/singularity) or conda

## documentation

Full documentation: [usage](./docs/usage.md).
Pipeline outputs: [output](./docs/output.md).

## acknowledgement

This pipeline is based on the Autocycler wrapper script by Ryan Wick: [autocycler_full.sh](https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh).

> [!NOTE]
> Please remember to properly cite Autocycler and the underlying tools it relies on when publishing results generated with this pipeline:
>
> Wick, R. R., Howden, B. P., & Stinear, T. P. (2025). Autocycler: long-read consensus assembly for bacterial genomes. Bioinformatics, 41(9). https://doi.org/10.1093/bioinformatics/btaf474
