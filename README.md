# F. pinicola annotation and comparative genomics pipeline

Snakemake workflow for repeat masking, gene prediction, functional annotation,
and three-species comparative genomic analysis based on a *Fomitopsis pinicola*
reference assembly.

This workflow was developed for the master's thesis *De novo assembly and
annotation of a near-chromosome-level reference genome for Fomitopsis pinicola,
and comparative genomic analysis with closely related species* (Nina Thorstensen,
University of South-Eastern Norway, 2026).

A companion workflow for assembly, polishing and heterozygosity analysis is
available 
[[LINK to assembly/QC repository here](https://github.com/ninkath/fpinicola-assembly-qc)]

## What this workflow does

Starting from a polished *F. pinicola* reference assembly (produced by the
companion workflow), the pipeline performs:

- **Repeat library construction** with RepeatModeler (LTR structural analysis
  enabled) combined with the Fungi subset of Dfam, used to soft-mask the genome
- **Gene prediction** with funannotate, integrating Augustus, SNAP, GlimmerHMM,
  GeneMark-ES and EVidenceModeler with external protein evidence
- **Functional annotation** with InterProScan, eggNOG-mapper (Basidiomycota
  scope), SignalP and dbCAN (HMMER + DIAMOND + dbCANsub consensus)
- **Comparative analysis** against two publicly available *Fomitopsis* genomes
  (*F. schrenkii* Fompi3 and *F. rosea* Fomro1, both from JGI MycoCosm),
  including BUSCO completeness, OrthoFinder orthology inference, COG functional
  category profiles, total and secreted CAZyme repertoires, pairwise macrosynteny and InterPro
  characterisation of species-specific orthogroups
- **Per-window tracks** used as additional input for the circos visualisation
  in the thesis (gene density, repeat composition by class)

A small number of supporting figures (circos plot, candidate MAT-B locus plot,
protein evidence file) are produced by manual scripts described below.

## Requirements

### Software

- **Snakemake** ≥ 9.13.4
- **Apptainer** (or Singularity) for containerised tools. The cache and tmp
  directories listed under `apptainer:` in `config.yaml`
  (default: `~/.apptainer/cache` and `~/.apptainer/tmp`) must exist before
  running the workflow.
- **Conda** or Micromamba for the conda environments under `envs/`:
  - `genemark_perl.yaml` provides the Perl modules required by GeneMark-ES
  - `plotting_python.yaml` provides Python with pandas and matplotlib
  - `r_plotting.yaml` provides R with ggplot2, pheatmap, ComplexUpset,
    and supporting packages

### Tools that must be installed locally

A few tools are not containerised and must be installed manually outside this
repository before running the workflow. The default paths in `config.yaml`
expect them under `tools/` and `secrets/`.

| Tool | Notes | Path in config |
|---|---|---|
| GeneMark-ES | Academic license required. The license key (`.gm_key`) must be placed at the path given in `annotation.genemark.genemark_key`. A separate Perl 5 directory with required modules is also expected. | `annotation.genemark.genemark_dir`, `annotation.genemark.genemark_key`, `annotation.genemark.perl5_dir` |
| InterProScan | Used for the reference annotation. The default expects v5.76-107.0 under `tools/interproscan-5.76-107.0/`. | `annotation.funannotate.interproscan_path` |
| SignalP | Academic license required. The default expects v5.0b under `tools/signalp-5.0b/`. | `annotation.signalp.dir`, `annotation.signalp.executable` |

### Reference databases

Several large reference databases must be downloaded and placed locally before
running the workflow. The default paths expect them under `resources/`. Note that some tools used have the capability to download standard databases when run. Information on how this is done can be found in the doumentation relevant to each program. This pipeline expects databases to be present before running.

| Database | Purpose | Path in config |
|---|---|---|
| funannotate database | Pfam, MEROPS, Augustus species models | `annotation.funannotate.db_dir` |
| eggNOG database | eggNOG-mapper functional annotation | `annotation.eggnog.db_dir` |
| dbCAN database | CAZyme HMM, DIAMOND and dbCANsub references | `annotation.dbcan.db_dir` |
| Dfam FamDB | Curated repeat families (Fungi subset extracted at runtime) | `annotation.dfam.db_dir` |
| BUSCO database | `polyporales_odb12` lineage dataset for the comparative BUSCO step | `qc.busco_db_path` |

The tool versions and database releases used for the thesis are listed in
Appendix E of the thesis.

### Protein evidence

Funannotate predict requires a protein evidence file. The script
`scripts/build_protein_evidence.sh` builds two evidence sets from local
input proteins:

- **Set A**: close relatives only (*F. betulina*, *F. rosea*, *Postia placenta*,
  *F. schrenkii* Fompi3)
- **Set B**: Set A plus the Swiss-Prot fungi subset (downloaded automatically)

Set B is the evidence file used for the thesis. The script expects raw input
proteomes under `proteins_raw/` and writes cleaned, ID-renamed evidence files
to `proteins_evidence/`. After the script has run, the active evidence file
must be placed at the path given in
`annotation.funannotate.protein_evidence.set_b` in `config.yaml`.

```bash
# Run from the repository root after placing input proteomes under proteins_raw/
bash scripts/build_protein_evidence.sh
```

### Outputs from the assembly/QC pipeline

A small number of files produced by the companion assembly/QC pipeline are
required as input for the circos plot and the candidate MAT-B locus plot
(see *Manual scripts* below). The expected location is
`external_inputs/assembly_pipeline/` in this repository.

## Repository contents

```
.
├── Snakefile                                  # Snakemake workflow definition
├── config.yaml                                # Sample, paths, parameters, container images
├── envs/
│   ├── genemark_perl.yaml                     # Perl modules for GeneMark-ES
│   ├── plotting_python.yaml                   # Python for table-generating scripts
│   └── r_plotting.yaml                        # R for plotting scripts
├── scripts/
│   ├── build_protein_evidence.sh              # Build protein evidence sets for funannotate
│   ├── parse_repeatmasker_out.py              # RepeatMasker .out → BED with repeat classes
│   ├── comparative_assembly_summary.py        # Combined assembly stats + BUSCO across species
│   ├── eggnog_cog_comparison.py               # COG category counts across species
│   ├── plot_cog_category_barplot.R            # COG barplot
│   ├── dbcan_class_summary.py                 # Per-class CAZyme counts across species
│   ├── dbcan_family_summary.py                # Per-family CAZyme counts across species
│   ├── plot_dbcan_family_heatmap.R            # CAZyme family heatmap
│   ├── merge_secreted_cazymes.py              # Intersect dbCAN consensus with SignalP
│   ├── secreted_dbcan_summary.py              # Per-family secreted CAZyme counts
│   ├── plot_secreted_dbcan_heatmap.R          # Secreted CAZyme heatmap
│   ├── orthofinder_presence_matrix.py         # Orthogroup presence/absence matrix
│   ├── orthogroup_intersection_summary.py     # Orthogroup set intersection counts
│   ├── orthogroup_species_specific_genes.py   # Genes from species-unique orthogroups
│   ├── orthogroup_summary_table.py            # Combined orthogroup summary table
│   ├── plot_orthogroup_upset.R                # Orthogroup UpSet plot
│   ├── species_specific_ipr_summary.py        # InterPro summary for species-specific orthogroups
│   ├── circos_genome_plot.R                   # Manual circos genome landscape figure
│   ├── mat_b_locus_plot.R                     # Manual candidate MAT-B locus figure
│   ├── synteny_extract_sequences.py           # Extract sequences for pairwise synteny
│   ├── synteny_prepare_plot.py                # NUCmer coordinate processing for synteny analysis
│   ├── synteny_plot_dotplot.R                 # Synteny figures
│   └── synteny_summary_table.py               # Summary statistics for the synteny analysis
├── README.md                                  # This file
└── LICENSE
```

The `results/`, `tools/`, `secrets/`, `resources/`, `proteins_raw/`,
`proteins_evidence/` and `external_inputs/` directories are not
version-controlled and are created or populated by the user.

## Configuration

The defaults in `config.yaml` are set to the values used in the thesis. Most
users will need to update at least the following:

```yaml
samples:
  fpindikaryon:                # Replace with your sample name

reference_sample: "fpindikaryon"

comparative:
  species:
    F_pinicola:
      proteins: "results/annotation/fpindikaryon/predict_results/fpindikaryon.proteins.fa"
      assembly: "results/assembly/fpindikaryon/assembly.fasta"
      busco_summary: "results/qc/busco_polished/fpindikaryon/short_summary.txt"

external:
  genome_windows_10kb: "external_inputs/assembly_pipeline/results/circos/fpindikaryon/genome_windows_10kb.bed"
```

The reference *F. pinicola* assembly produced by the assembly/QC pipeline must
be placed at the path given by `comparative.species.F_pinicola.assembly`
(`results/assembly/{reference_sample}/assembly.fasta` by default).
In the thesis project, this was the polished Medaka consensus copied in and
renamed to `assembly.fasta`.

## Running the workflow

This workflow does **not** define a default `rule all`. Each result is generated
by specifying its target output file. Rules are intended to be run in the order
described in the thesis methods.

A typical run sequence:

```bash
# 1. Repeat library construction and masking
snakemake --use-apptainer --cores 16 \
    results/annotation/fpindikaryon/masked_genome.fasta

# 2. Ab initio gene prediction (GeneMark-ES) and consensus prediction (funannotate)
snakemake --use-apptainer --use-conda --cores 16 \
    results/annotation/fpindikaryon/predict_results/fpindikaryon.proteins.fa

# 3. Functional annotation
snakemake --use-apptainer --cores 8 \
    results/annotation/fpindikaryon/eggnog/annotations.emapper.annotations \
    results/annotation/fpindikaryon/interproscan/interproscan.tsv \
    results/annotation/fpindikaryon/signalp/fpindikaryon_signalp_summary.txt \
    results/annotation/fpindikaryon/dbcan/overview.tsv

# 4. Final funannotate annotation set
snakemake --use-apptainer --cores 8 \
    results/annotation/fpindikaryon/final_annotation.done

# 5. Comparative analysis (run from this point only after the comparative
#    species protein and assembly files are placed under comparative/)
snakemake --use-apptainer --use-conda --cores 16 \
    results/comparative/tables/comparative_assembly_summary.tsv \
    results/comparative/tables/orthogroup_summary_table.tsv \
    results/comparative/tables/orthogroup_intersection_summary.tsv \
    results/comparative/tables/species_specific_orthogroups.tsv \
    results/comparative/tables/cog_category_comparison.tsv \
    results/comparative/tables/dbcan_class_counts.tsv \
    results/comparative/tables/dbcan_family_counts.tsv \
    results/comparative/tables/secreted_dbcan_family_counts.tsv \
    results/comparative/tables/species_specific_ipr_summary.tsv

# 6. Comparative plots
snakemake --use-apptainer --use-conda --cores 4 \
    results/comparative/plots/orthogroups_upset.png \
    results/comparative/plots/cog_category_barplot.png \
    results/comparative/plots/dbcan_family_heatmap.png \
    results/comparative/plots/secreted_dbcan_family_heatmap.png

snakemake --use-apptainer --use-conda --cores 4 synteny_dotplot_all

# 7. Per-window tracks for the circos plot
snakemake --use-apptainer --cores 4 \
    results/circos/fpindikaryon/gene_density_10kb.tsv \
    results/circos/fpindikaryon/repeat_bp_density_10kb.tsv
```

A dry-run of any command above (add `-n`) shows the rules that will execute
without producing files.

## Manual scripts

A few scripts are run manually outside Snakemake. All of them use relative
paths and assume the repository root as the working directory. Before running,
open the script and confirm that the paths and the `REFERENCE` variable match
your setup.

### Building the protein evidence file

```bash
# From the repository root, after placing raw proteomes under proteins_raw/
bash scripts/build_protein_evidence.sh
```

### Circos genome landscape figure

The circos plot integrates outputs from both this pipeline (gene density,
repeat tracks) and from the assembly/QC pipeline (GC, coverage, telomere
windows, heterozygosity windows). The companion pipeline outputs are expected
under `external_inputs/assembly_pipeline/`. Check input file paths before running.

```bash
# From the repository root, after the dependencies above are in place
Rscript scripts/circos_genome_plot.R
```

### Candidate MAT-B locus figure

This script reads the funannotate GFF for the reference and the per-window
heterozygosity track from the assembly/QC pipeline.

```bash
Rscript scripts/mat_b_locus_plot.R
```

## License

MIT License (see LICENSE).
