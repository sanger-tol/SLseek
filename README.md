# SL discovery

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)

🧱 Add links to singularity and docker

## Introduction
**SLseek** is a bioinformatics pipeline for the discovery of highly represented sequences at the ends of transcripts, usually spliced leader (SL) sequences in organisms with trans-splicing. It analyses full-length transcriptomic data using a k-mer-based approach, avoiding costly steps like mapping or transcriptome assembly, and does not require a reference genome. 

This pipeline is built using [Nextflow](https://www.nextflow.io), a scalable workflow management tool. 

## Pipeline summary

1. Extract ends of transcripts from the input full-length transcripts;
2. Count kmers using either [FastK](https://github.com/thegenemyers/FASTK) and/or [Jellyfish](https://github.com/gmarcais/Jellyfish);
3. Get Putative SL sequences (e.g, sequences composed of kmers highly present at the ends of transcripts, within a specified size and mostly non-repetitive);
4. Cluster Putative SLs using [CD-HIT](https://github.com/weizhongli/cdhit).

Optionally, the pipeline can generate a k-mer histogram and a heatmap plot with putative SL assessment metrics. 

![](https://github.com/Beatriz-Estevam/SLseek/blob/main/figures/SLseek_representation.png)

## Usage

> [!NOTE]
> Please refer to [this page](https://nf-co.re/docs/usage/installation) to install and set-up Nextflow. You can test the installation following [this link](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test`.

🧱 Add links to singularity and docker

You can run the pipeline using:

```bash
nextflow
```





