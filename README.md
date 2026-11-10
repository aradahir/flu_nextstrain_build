# Nextstrain build for human influenza surveilance

Phylogenetic analysis of influenza sequences using Nextstrain tools ([Augur](https://docs.nextstrain.org/projects/augur) and [Auspice](https://docs.nextstrain.org/projects/auspice)). The workflow is adapted from the Nextstrain [seasonal-flu](https://github.com/nextstrain/seasonal-flu) repository.

This README covers two ways to run the build:

1. [Step-by-step Augur commands](#step-by-step-build) (easiest to understand and debug)
2. [Snakemake workflow](#snakemake-workflow) (to run at once)

## Contents

- [Installation](#installation)
- [Repository layout](#repository-layout)
- [Input data](#input-data)
- [Step-by-step build](#step-by-step-build)
- [Snakemake workflow](#snakemake-workflow)
- [Optional additions](#optional-additions)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements](#acknowledgements)

## Installation

Install the Nextstrain CLI and the tools it needs by following the [Nextstrain installation guide](https://docs.nextstrain.org/en/latest/install.html).

For a conda setup:

```bash
conda env create -f environment.yml
conda activate nextstrain_build
```

Check the install:

```bash
augur --version
nextstrain --version
```

## Repository layout

```
.
├── snakefile                 # Snakemake workflow
├── nexclade/                 # Automatic clade assignment and haplotype source
├── config/
│   ├── metadata/             # reference backbone epidemiology information
│   ├── reference/            # reference database (Genbank)
│   ├── reference_fasta/      # reference fasta
│   ├── vaccine.json          # vaccine strain to highlight (must be in the pattern `{subtype}_vaccine.json`)
│   ├── colors.tsv            # custom colours for Auspice
│   └── auspice_config.json   # Auspice display settings
├── data/
│   ├── concatenated_files    # input sequences with reference sequences (created by the build)
│   ├── sequences.fasta       # input sequences (must be in the pattern `{subtype}_{segment}.fasta`)
│   └── metadata.tsv          # input metadata (must be in the pattern `{subtype}_{segment}_metadata.tsv`)
├── script/                   # helper scripts (frequencies, haplotypes)
├── logs/                     # log files for debugging
├── results/                  # intermediate files (created by the build)
└── auspice/                  # final JSONs for Auspice (created by the build)
```

## Input data

| File | Description |
|---|---|
| `data/sequences.fasta` | HA (or NA segment) nucleotide sequences. |
| `data/metadata.tsv` | Tab-separated metadata. Must include `strain` and `date`. The `strain` values must match the FASTA headers exactly. Other columns used here: `country`, `region`, `clade`, `subclade`, `location`, `passage`, `state`. |

## Step-by-step build

These steps were used for creating the pipeline from the repository.

```mermaid
flowchart LR
    A[Concatenate<br/> with reference files] --> B["Clade calling<br/>(nextclade)"]
    B --> C[Assign haplotype]
    C --> D[Index]
    D --> E["Align"]
    E --> F["Tree"]
    F --> G["Refine Tree<br/>time tree"]
    G --> H["Traits"]
    H --> I["Ancestral"]
    I --> J["Translate"]
    J --> K["Estimate Frequency"]
    K --> L["Export"]
    L --> M[View]


```

Then open the address printed in the terminal (usually http://localhost:4000). Alternatively, drag `auspice/{subtype}_{segment}.json` onto [auspice.us](https://auspice.us).

## Snakemake workflow

The `Snakefile` runs the same steps for every `{subtype}_{segment}` combination, using wildcards. For example, to build h1n1pdm_HA:

```bash
snakemake --cores 4 
```

Preview what would run without executing it:

```bash
snakemake -n --cores 4 
```

View the result:

```bash
nextstrain view results/h1n1pdm_ha/auspice/
```

## Optional additions

**Tip frequencies.** `script/estimate_freq.py` estimates subclade frequencies with a kernel density estimate and writes `*_tip-frequencies.json` next to the Auspice JSON. Auspice loads it automatically when it has the same prefix as the main JSON. Do not pass it to `--node-data`.

**Derived haplotypes.** `script/make_haplotype.py` adds a `haplotype` column to the metadata, built from amino acid substitutions that go beyond each subclade's defining mutations. Add `haplotype` to the `filters` and `color_options` in `config/auspice_config.json`.

**Vaccine strains.** To mark vaccine strains with an X on the tree, create a node-data JSON (for example `config/{subtype}_vaccine.json`) with a `vaccine` attribute for each strain name, and add it to `--node-data` in the export step.

## Troubleshooting

| Problem | Likely cause |
|---|---|
| `Illegal trailing comma before end of object` | A trailing comma in `config/auspice_config.json`. Check with `python -m json.tool config/auspice_config.json`. |
| `has empty or nonexistent nodes and branches` | A tip-frequencies JSON was passed to `--node-data`. Remove it from that list. |
| Strains missing from the tree | `strain` names in the metadata do not match the FASTA headers. |
| All amino acid mutations empty | The gene names in the reference GenBank do not match what `augur translate` expects. |

## Acknowledgements

This build is adapted from the [Nextstrain seasonal-flu workflow](https://github.com/nextstrain/seasonal-flu). Please also cite Augur and Auspice:

> Hadfield J, Megill C, Bell SM, et al. Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics* 2018.