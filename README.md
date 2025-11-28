# SL discovery

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## 🚀 Introduction
**SLseek** is a bioinformatics pipeline for discovering highly represented sequences at the ends of transcripts, usually spliced leader (SL) sequences in organisms with trans-splicing. It analyses full-length transcriptomic data using a k-mer-based approach, avoiding costly steps such as mapping or transcriptome assembly, and does not require a reference genome. 

This pipeline is built using [Nextflow](https://www.nextflow.io), a scalable workflow management tool. 

## 🔁 Pipeline summary

1. Extract ends of transcripts from the input full-length transcripts;
2. Count kmers using either [FastK](https://github.com/thegenemyers/FASTK) and/or [Jellyfish](https://github.com/gmarcais/Jellyfish);
3. Get Putative SL sequences (e.g, sequences composed of kmers highly present at the ends of transcripts, within a specified size and mostly non-repetitive);
4. Cluster Putative SLs using [CD-HIT](https://github.com/weizhongli/cdhit).

Optionally, the pipeline can generate a k-mer histogram and a heatmap plot with putative SL assessment metrics. 

![](https://github.com/Beatriz-Estevam/SLseek/blob/main/figures/SLseek_representation.png)

## ⚙️ Usage

> [!NOTE]
> Please refer to [this page](https://nf-co.re/docs/usage/installation) to install and set-up Nextflow. You can test the installation following [this link](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test`.

```bash
nextflow run transsplicing.nf -profile [standard | lsf [ ,singularity ]] \
    --outdir <path/to/outdir/> \
    --transcripts <path/to/transcripts.fasta> \
    --kmer_counting_tool (fastk | jellyfish) \
    -c [transsplicing.config] \
    [ -resume ]
```

## 📖 Understanding parameters 
| Parameter | Default | Required | Description | Options/Input Type |
| :--- | :--- | :--- | :--- | :--- |
| `outdir` | | **Yes** | A path to where saving outputs. | `string` |
| `transcripts` |  | **Yes** | A fasta file with full-length transcripts | `string` |
| `length` | `100` | **No**  | Length of extracted sequences for analysis - Should encompass SL size. | `int` |
| `k` | `21` | **No** | K-mer size for analysis - Should be less than SL-size and typically an odd number. | `int` |
| `lower_cov_thrsld_counting` | `150` | **No** | Minimum coverage of kmers to filter using Jellyfish or FastK | `int` |
| `lower_cov_thrsld_extracting` | `150` | **No** | Minimum coverage of kmers to filter in k-mer analysis steps | `int` |
| `max_distance` | `3`  | **No** | Minimum distance to consider a K-mer within the same group | `int` |
| `extra_border` | `0` | **No** | Number of extranucleotides flanking putative SLs from K-mers. | `int` |
| `size_limit_max` | `30` | **No** | Maximum size of a Putative SL. | `int` |
| `size_limit_min` | `20` | **No** | Miiumum size of a Putative SL. | `int` |
| `entropy_lim` | `1.40` | **No** | Minumum entropy value for Putative SL to avoid repetitive sequences | `float` |
| `nucleotide_limit` | **0.7** | **No** | Percentage of nucleotides to consider a kmer as uninformative. E.g., for 0.7, kmers with 70% of the same nucleotides will be removed | `float` |
| `identity_thrsld` | `0.90` | **No** | Identity threshold for CD-HIT Putative SL clustering. | `float` |
| `heatmap` | `false` | **No** | If called or true, will generate a heatmap plot with putative SL analysis metrics | `boolean` (E.g., true; false) |
| `kmerplot` | `false` | **No** | If called or true, will generate a kmer histogram with kmer counting process | `boolean` (E.g., true; false) |
| `kmer_counting_tool` | `jellyfish`  | **Yes** | kmer counting tool | `string` (E.g., jellyfish or fastk)|


## 🧪 Test data

A  set of data to test is provided within the project repository to help users quickly verify the installation and functionality of the pipeline. You can find all relevant files in the dedicated `test/` directory of this repository. 

### Running tests 

! Assuming you are in the root directory of the pipeline


```bash
nextflow run transsplicing.nf \
        -c transsplicing.config
        --transcripts "$INPUT_PATH" \
        --k 17 \
        --outdir "test/results/" \
        --kmer_counting_tool "fastk" \
        --heatmap --kmerplot \
        -profile lsf, singularity \
        -resume 
```

**RESULTS**
- Kmer plot
![Descrição da imagem](test/results//example.png)

- Heatmap


- Main result

```bash
nextflow run transsplicing.nf \
        -c transsplicing.config
        --transcripts "$INPUT_PATH" \
        --lower_cov_thrsld_counting 500 \
        --k 17 \
        --outdir "test/results/" \
        --kmer_counting_tool "jellyfish" \
        --heatmap --kmerplot \
        -profile lsf, singularity \
        -resume 
```


```bash
nextflow run transsplicing.nf \
        -c transsplicing.config
        --transcripts "$INPUT_PATH" \
        --lower_cov_thrsld_counting 500 \
        --k 21 \
        --outdir "test/results/" \
        --kmer_counting_tool "fastk" \
        --heatmap --kmerplot \
        -profile lsf, singularity \
        -resume 
```

**RESULTS**
- Kmer plot

- Heatmap


- Main result




