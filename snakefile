import os
from pathlib import Path

configfile: "config.yaml"

REFERENCE = config["reference_sample"]
COMPARATIVE = config["comparative"]["species"]

# GeneMark paths
GENEMARK_DIR = os.path.expanduser(config["annotation"]["genemark"]["genemark_dir"])
GENEMARK_KEY = os.path.expanduser(config["annotation"]["genemark"]["genemark_key"])

# No default rule_all target is defined in this workflow.
# Run individual rules manually in the desired order.

# ============================================================================
# HELPER FUNCTIONS FOR COMPARATIVE GENOMICS
# ============================================================================
def dbcan_overview_for_species(wc):
    if wc.species == "F_pinicola":
        return f"results/annotation/{REFERENCE}/dbcan/overview.tsv"
    return f"results/comparative/dbcan/{wc.species}/overview.tsv"

def signalp_summary_for_species(wc):
    if wc.species == "F_pinicola":
        return f"results/annotation/{REFERENCE}/signalp/{REFERENCE}_signalp_summary.txt"
    return f"results/comparative/signalp/{wc.species}/{wc.species}_signalp_summary.txt"

def eggnog_for_species(wc):
    if wc.species == "F_pinicola":
        return f"results/annotation/{REFERENCE}/eggnog/annotations.emapper.annotations"
    return f"results/comparative/eggnog/{wc.species}/annotations.emapper.annotations"

def interproscan_tsv_for_species(wc):
    if wc.species == "F_pinicola":
        return f"results/annotation/{REFERENCE}/interproscan/interproscan.tsv"
    return COMPARATIVE[wc.species]["ipr"]

def busco_summary_for_species(wc):
    return COMPARATIVE[wc.species]["busco_summary"]

# ============================================================================
# ANNOTATION (on reference genome)
# ============================================================================
# repeatmodeler_builddb: build BLAST database from the assembly for RepeatModeler
rule repeatmodeler_builddb:
    input:
        assembly = f"results/assembly/{REFERENCE}/assembly.fasta"
    output:
        nhr = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nhr",
        nin = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nin",
        nsq = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nsq"
    params:
        outdir = f"results/annotation/{REFERENCE}/repeatmodeler",
        dbname = "sample_db",
        cache = os.path.expanduser(config["apptainer"]["cache_dir"]),
        tmpdir = os.path.expanduser(config["apptainer"]["tmp_dir"]),
        image = config["containers"]["repeatmodeler"]
    threads: config["threads"]["light"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}

        export APPTAINER_CACHEDIR="{params.cache}"
        export APPTAINER_TMPDIR="{params.tmpdir}"

        INPUT_ASM=$(realpath {input.assembly})
        OUTDIR=$(realpath {params.outdir})

        cd "$OUTDIR"

        apptainer exec --bind "$PWD" "{params.image}" \
            BuildDatabase -name {params.dbname} "$INPUT_ASM"

        test -s {params.dbname}.nhr
        test -s {params.dbname}.nin
        test -s {params.dbname}.nsq
        """

# repeatmodeler_run: de novo repeat library construction (LTR structural analysis enabled)
rule repeatmodeler_run:
    input:
        nhr = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nhr",
        nin = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nin",
        nsq = f"results/annotation/{REFERENCE}/repeatmodeler/sample_db.nsq"
    output:
        lib = f"results/annotation/{REFERENCE}/repeatmodeler/consensi.fa.classified"
    params:
        outdir = f"results/annotation/{REFERENCE}/repeatmodeler",
        dbname = "sample_db",
        srand=config["annotation"]["repeatmodeler"]["srand"],
        enable_ltr = config["annotation"]["repeatmodeler"]["enable_ltr"],
        cache = os.path.expanduser(config["apptainer"]["cache_dir"]),
        tmpdir = os.path.expanduser(config["apptainer"]["tmp_dir"]),
        image = config["containers"]["repeatmodeler"]
    threads: config["threads"]["heavy"]
    shell:
        r"""
        set -euo pipefail

        export APPTAINER_CACHEDIR="{params.cache}"
        export APPTAINER_TMPDIR="{params.tmpdir}"

        OUTDIR=$(realpath {params.outdir})
        cd "$OUTDIR"

        if [ ! -s "{params.dbname}-families.fa" ]; then
            find . -maxdepth 1 -type d -name "RM_*" -exec rm -rf {{}} + || true

            if [ "{params.enable_ltr}" = "True" ]; then
                echo "Running RepeatModeler WITH LTR analysis"
                apptainer exec --bind "$PWD" "{params.image}" \
                    RepeatModeler -database {params.dbname} -threads {threads} -LTRStruct -srand {params.srand} \
                    || {{
                        echo "LTR analysis failed, retrying without LTRStruct"
                        apptainer exec --bind "$PWD" "{params.image}" \
                            RepeatModeler -database {params.dbname} -threads {threads}
                    }}
            else
                echo "Running RepeatModeler WITHOUT LTR analysis"
                apptainer exec --bind "$PWD" "{params.image}" \
                    RepeatModeler -database {params.dbname} -threads {threads}
            fi
        fi

        if [ -s "{params.dbname}-families.fa" ]; then
            cp "{params.dbname}-families.fa" "consensi.fa.classified"
        else
            echo "ERROR: No consensus sequences found"
            exit 1
        fi

        test -s "consensi.fa.classified"
        echo "RepeatModeler complete"
        """

# extract_dfam_fungi: extract Fungi-specific repeats from local Dfam FamDB
rule extract_dfam_fungi:
    input:
        db_dir = config["annotation"]["dfam"]["db_dir"]
    output:
        fasta = f"results/annotation/{REFERENCE}/repeatmodeler/dfam_fungi.fasta"
    container: config["containers"]["tetools"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})
        
        /opt/RepeatMasker/famdb.py -i {input.db_dir} families -f fasta_name -ad Fungi > {output.fasta}
        
        test -s {output.fasta}
        """

# combine_repeat_libraries: concatenate Dfam and de novo libraries into a single FASTA
rule combine_repeat_libraries:
    input:
        dfam = f"results/annotation/{REFERENCE}/repeatmodeler/dfam_fungi.fasta",
        denovo = f"results/annotation/{REFERENCE}/repeatmodeler/consensi.fa.classified"
    output:
        combined = f"results/annotation/{REFERENCE}/repeatmodeler/combined_repeats.fasta"
    shell:
        r"""
        set -euo pipefail
        
        cat {input.dfam} {input.denovo} > {output.combined}
        test -s {output.combined}
        """

# repeatmasker_mask: soft-mask the assembly using the combined repeat library
rule repeatmasker_mask:
    input:
        assembly = f"results/assembly/{REFERENCE}/assembly.fasta",
        repeats = f"results/annotation/{REFERENCE}/repeatmodeler/combined_repeats.fasta"
    output:
        masked = f"results/annotation/{REFERENCE}/masked_genome.fasta"
    params:
        outdir = f"results/annotation/{REFERENCE}",
        use_repeatmodeler = config["annotation"]["repeatmodeler"]["use_with_funannotate"],
        cache = os.path.expanduser(config["apptainer"]["cache_dir"]),
        tmpdir = os.path.expanduser(config["apptainer"]["tmp_dir"]),
        rm_image = config["containers"]["tetools"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}

        export APPTAINER_CACHEDIR="{params.cache}"
        export APPTAINER_TMPDIR="{params.tmpdir}"

        ASSEMBLY=$(realpath {input.assembly})
        REPEATLIB=$(realpath {input.repeats})
        OUTDIR=$(realpath {params.outdir})

        cd "$OUTDIR"

        if [ "{params.use_repeatmodeler}" = "True" ]; then
            echo "Masking with RepeatMasker using combined library (Dfam + de novo)..."
            apptainer exec --bind "$PWD" "{params.rm_image}" \
                RepeatMasker -lib "$REPEATLIB" \
                    -pa {threads} \
                    -xsmall \
                    -dir . \
                    "$ASSEMBLY"

            BASENAME=$(basename "$ASSEMBLY")
            mv "$BASENAME.masked" "masked_genome.fasta"
        else
            echo "Copying assembly without masking..."
            cp "$ASSEMBLY" "masked_genome.fasta"
        fi
        """

# genemark_es: ab initio gene prediction with GeneMark-ES in fungal mode
rule genemark_es:
    input:
        masked = f"results/annotation/{REFERENCE}/masked_genome.fasta",
        gm_key = GENEMARK_KEY
    output:
        gtf = f"results/annotation/{REFERENCE}/genemark/genemark.gtf",
        done = f"results/annotation/{REFERENCE}/genemark/genemark.done"
    params:
        outdir = f"results/annotation/{REFERENCE}/genemark",
        gm_dir = GENEMARK_DIR,
        perl5_dir = config["annotation"]["genemark"]["perl5_dir"]
    threads: config["threads"]["medium"]
    conda: "envs/genemark_perl.yaml"
    shell:
        r"""
        set -euo pipefail

        GM_DIR="$(realpath {params.gm_dir})"
        GENEMARK_SCRIPT="$GM_DIR/gmes_petap.pl"
        INPUT_FASTA="$(realpath {input.masked})"
        PERL5_DIR="$(realpath {params.perl5_dir})"

        export PERL5LIB="$PERL5_DIR/lib/perl5"

        export TMPDIR="{params.outdir}/tmp"
        mkdir -p "$TMPDIR"

        mkdir -p "$HOME"
        ln -sf "$(realpath {input.gm_key})" "$HOME/.gm_key"

        export PATH="$GM_DIR:$GM_DIR/ProtHint/bin:$PATH"

        perl -MHash::Merge -e 'print "Hash::Merge found\n"' || {{
            echo "ERROR: Hash::Merge module not found"
            exit 1
        }}

        mkdir -p {params.outdir}
        cd {params.outdir}

        perl "$GENEMARK_SCRIPT" \
            --ES --fungus \
            --cores {threads} \
            --max_intron 3000 \
            --soft_mask 2000 \
            --sequence "$INPUT_FASTA"

        if [ -f "genemark.gtf" ]; then
            :
        elif [ -f "output/gene.gtf" ]; then
            cp output/gene.gtf genemark.gtf
        else
            echo "ERROR: Could not find GeneMark GTF output" >&2
            ls -R
            exit 1
        fi

        test -s "genemark.gtf"
        touch genemark.done
        """

# funannotate_predict: consensus gene prediction (Augustus, SNAP, GlimmerHMM, GeneMark, EVidenceModeler)
rule funannotate_predict:
    input:
        masked = f"results/annotation/{REFERENCE}/masked_genome.fasta",
        genemark_done = f"results/annotation/{REFERENCE}/genemark/genemark.done",
        genemark_gtf = f"results/annotation/{REFERENCE}/genemark/genemark.gtf"
    output:
        done = f"results/annotation/{REFERENCE}/predict.done",
        gbk = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.gbk",
        proteins = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.proteins.fa",
        gff = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.gff3"
    params:
        outdir = f"results/annotation/{REFERENCE}",
        species = config["annotation"]["species_name"],
        augustus_species = config["annotation"]["augustus_species"],
        busco_seed_species = config["annotation"].get("busco_seed_species", ""),
        busco_db = config["annotation"]["funannotate"]["busco_db"],
        db_dir = config["annotation"]["funannotate"]["db_dir"],
        isolate = REFERENCE,
        protein_evidence_path = (
            config["annotation"]["funannotate"]["protein_evidence"].get(
                config["annotation"]["funannotate"]["protein_evidence"].get("active", ""),
                ""
            ) if config["annotation"]["funannotate"].get("protein_evidence") else ""
        ),
        weights = "augustus:1 HiQ:2 GeneMark:1 snap:1 glimmerhmm:1 proteins:1"
    threads: config["threads"]["heavy"]
    log: f"logs/funannotate/{REFERENCE}_predict.log"
    container: config["containers"]["funannotate"]
    shell:
        r"""
        set -euo pipefail

        MASKED=$(realpath {input.masked})
        GENEMARK_GTF=$(realpath {input.genemark_gtf})
        OUTDIR=$(realpath {params.outdir})
        DB_DIR=$(realpath {params.db_dir})

        mkdir -p "$OUTDIR"
        mkdir -p "$(dirname {log})"

        export FUNANNOTATE_DB="$DB_DIR"

        echo "=== funannotate predict ===" | tee {log}
        echo "Species: {params.species}" | tee -a {log}
        echo "Augustus species: {params.augustus_species}" | tee -a {log}
        echo "BUSCO seed species: {params.busco_seed_species}" | tee -a {log}
        echo "BUSCO db: {params.busco_db}" | tee -a {log}
        echo "EVM weights: {params.weights}" | tee -a {log}

        PROT_ARG=""
        if [ -n "{params.protein_evidence_path}" ]; then
            if [ -f "{params.protein_evidence_path}" ]; then
                PROT_FILE=$(realpath "{params.protein_evidence_path}")
                PROT_ARG="--protein_evidence $PROT_FILE"
                echo "Using protein evidence: $PROT_FILE" | tee -a {log}
            else
                echo "WARNING: Protein evidence file not found: {params.protein_evidence_path}" | tee -a {log}
                echo "Continuing without protein evidence..." | tee -a {log}
            fi
        else
            echo "No protein evidence configured" | tee -a {log}
        fi

        BUSCO_SEED_ARG=""
        if [ -n "{params.busco_seed_species}" ]; then
            BUSCO_SEED_ARG="--busco_seed_species {params.busco_seed_species}"
            echo "Using BUSCO seed: {params.busco_seed_species}" | tee -a {log}
        fi

        funannotate predict \
            -i "$MASKED" \
            -o "$OUTDIR" \
            --species "{params.species}" \
            --isolate "{params.isolate}" \
            --augustus_species {params.augustus_species} \
            $BUSCO_SEED_ARG \
            --busco_db {params.busco_db} \
            --genemark_gtf "$GENEMARK_GTF" \
            --weights {params.weights} \
            $PROT_ARG \
            --ploidy 1 \
            --cpus {threads} \
            2>&1 | tee -a {log}

        PRED_DIR="$OUTDIR/predict_results"

        link_one () {{
          local pattern="$1"
          local dest="$2"
          local src
          src="$(ls -1 "$PRED_DIR"/$pattern 2>/dev/null | sort | head -n1 || true)"
          if [ -n "${{src:-}}" ] && [ -s "$src" ]; then
            ln -sfn "$(basename "$src")" "$dest" || cp -f "$src" "$dest"
          else
            echo "ERROR: Could not find file matching '$pattern' in $PRED_DIR" >&2
            ls -1 "$PRED_DIR" || true
            exit 1
          fi
        }}

        link_one "*.gbk" {output.gbk}
        link_one "*proteins.fa" {output.proteins}
        link_one "*.gff3" {output.gff}

        touch {output.done}
        """

# eggnog_map: functional annotation with eggNOG-mapper, taxonomic scope Basidiomycota
rule eggnog_map:
    input:
        faa = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.proteins.fa"
    output:
        ann = f"results/annotation/{REFERENCE}/eggnog/annotations.emapper.annotations"
    params:
        data_dir = config["annotation"]["eggnog"]["db_dir"],
        outdir = f"results/annotation/{REFERENCE}/eggnog",
        prefix = "annotations",
        tax = config["annotation"]["eggnog"]["tax_scope"]
    threads: config["threads"]["medium"]
    container: config["containers"]["eggnog-mapper"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        emapper.py \
            -i {input.faa} \
            --itype proteins \
            -o {params.prefix} \
            --output_dir {params.outdir} \
            --cpu {threads} \
            --usemem \
            --tax_scope {params.tax} \
            --data_dir "{params.data_dir}"
        """

# interproscan: protein domain and family annotation with InterProScan
rule interproscan:
    input:
        proteins = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.proteins.fa"
    output:
        tsv = f"results/annotation/{REFERENCE}/interproscan/interproscan.tsv",
        xml = f"results/annotation/{REFERENCE}/interproscan/interproscan.xml"
    params:
        iprscan_sh = config["annotation"]["funannotate"]["interproscan_path"],
        outdir = f"results/annotation/{REFERENCE}/interproscan",
        tmpdir = f"results/annotation/{REFERENCE}/interproscan/temp"
    threads: config["threads"]["medium"]
    resources:
        mem_mb = 16000
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}

        export JAVA_OPTS="-Xms2g -Xmx14g"
        export TMPDIR={params.tmpdir}

        if [ ! -s {input.proteins} ]; then
            echo "ERROR: Input file is empty or missing: {input.proteins}"
            exit 1
        fi

        {params.iprscan_sh} \
            -i {input.proteins} \
            -f tsv,xml \
            -cpu 4 \
            -b {params.outdir}/interproscan \
            -T {params.tmpdir} \
            -dp \
            -goterms

        test -f {output.tsv}
        test -f {output.xml}

        rm -rf {params.tmpdir} 
        """

# signalp_reference: signal peptide prediction for the reference proteins
rule signalp_reference:
    input:
        proteins = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.proteins.fa"
    output:
        summary = f"results/annotation/{REFERENCE}/signalp/{REFERENCE}_signalp_summary.txt",
        done = f"results/annotation/{REFERENCE}/signalp/signalp.done"
    params:
        outdir = f"results/annotation/{REFERENCE}/signalp",
        signalp_dir = config["annotation"]["signalp"]["dir"],
        organism = "euk"
    threads: 4
    log:
        f"logs/signalp/{REFERENCE}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname {log})"

        OUTDIR=$(realpath "{params.outdir}")
        PROTEINS=$(realpath "{input.proteins}")
        SIGNALP_DIR=$(realpath "{params.signalp_dir}")
        LOGFILE=$(realpath "{log}")
        SUMMARY_OUT=$(realpath "$(dirname {output.summary})")/$(basename {output.summary})
        DONE_OUT=$(realpath "$(dirname {output.done})")/$(basename {output.done})

        cd "$SIGNALP_DIR/bin"

        ./signalp \
            -fasta "$PROTEINS" \
            -org {params.organism} \
            -format short \
            -prefix "$OUTDIR/{REFERENCE}_signalp" \
            2>&1 | tee "$LOGFILE"

        if [ -f "$OUTDIR/{REFERENCE}_signalp_summary.signalp5" ]; then
            cp "$OUTDIR/{REFERENCE}_signalp_summary.signalp5" "$SUMMARY_OUT"
        else
            echo "ERROR: Expected SignalP summary file not found in $OUTDIR" >&2
            ls -lah "$OUTDIR" >&2 || true
            exit 1
        fi

        test -s "$SUMMARY_OUT"
        touch "$DONE_OUT"
        """

# dbcan_reference: CAZyme annotation with dbCAN (HMMER, DIAMOND, dbCANsub)
rule dbcan_reference:
    input:
        proteins = f"results/annotation/{REFERENCE}/predict_results/{REFERENCE}.proteins.fa"
    output:
        overview = f"results/annotation/{REFERENCE}/dbcan/overview.tsv",
        hmm = f"results/annotation/{REFERENCE}/dbcan/dbCAN_hmm_results.tsv",
        sub = f"results/annotation/{REFERENCE}/dbcan/dbCANsub_hmm_results.tsv",
        diamond = f"results/annotation/{REFERENCE}/dbcan/diamond.out",
        done = f"results/annotation/{REFERENCE}/dbcan/dbcan.done"
    params:
        outdir = f"results/annotation/{REFERENCE}/dbcan",
        db_dir = config["annotation"]["dbcan"]["db_dir"],
        evalue = config["annotation"]["dbcan"]["evalue"],
        coverage = config["annotation"]["dbcan"]["coverage"]
    threads: config["threads"]["medium"]
    log:
        f"logs/dbcan/{REFERENCE}.log"
    container:
        config["containers"]["dbcan"]
    shell:
        r"""
        set -euo pipefail

        PROTEINS=$(realpath {input.proteins})
        mkdir -p {params.outdir}
        mkdir -p "$(dirname {log})"

        OUTDIR=$(realpath {params.outdir})
        DB_DIR=$(realpath {params.db_dir})

        run_dbcan CAZyme_annotation \
            --input_raw_data "$PROTEINS" \
            --output_dir "$OUTDIR" \
            --db_dir "$DB_DIR" \
            --mode protein \
            --methods hmm --methods diamond --methods dbCANsub \
            --e_value_threshold_dbcan {params.evalue} \
            --coverage_threshold_dbcan {params.coverage} \
            --e_value_threshold_dbsub {params.evalue} \
            --coverage_threshold_dbsub {params.coverage} \
            --threads {threads} \
            2>&1 | tee {log}

        test -f {output.overview}
        touch {output.done}
        """

# funannotate_annotate: combine all annotation outputs into the final annotated gene set
rule funannotate_annotate:
    input:
        predict = f"results/annotation/{REFERENCE}/predict.done",
        interproscan = f"results/annotation/{REFERENCE}/interproscan/interproscan.xml",
        eggnog = f"results/annotation/{REFERENCE}/eggnog/annotations.emapper.annotations",
        signalp = f"results/annotation/{REFERENCE}/signalp/{REFERENCE}_signalp_summary.txt",
        dbcan = f"results/annotation/{REFERENCE}/dbcan/overview.tsv"
    output:
        done = f"results/annotation/{REFERENCE}/final_annotation.done",
        gff = f"results/annotation/{REFERENCE}/annotate_results/{REFERENCE}.gff3",
        gbk = f"results/annotation/{REFERENCE}/annotate_results/{REFERENCE}.gbk",
        proteins = f"results/annotation/{REFERENCE}/annotate_results/{REFERENCE}.proteins.fa"
    params:
        outdir = f"results/annotation/{REFERENCE}",
        db_dir = config["annotation"]["funannotate"]["db_dir"],
        busco_db = config["annotation"]["funannotate"]["busco_db"],
        container = config["containers"]["funannotate"],
        ref_name = REFERENCE
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail

        DB_PATH=$(realpath {params.db_dir})
        OUT_DIR=$(realpath {params.outdir})
        EGGNOG=$(realpath {input.eggnog})
        IPRSCAN=$(realpath {input.interproscan})
        SIGNALP=$(realpath {input.signalp})
        DBCAN=$(realpath {input.dbcan})

        mkdir -p "$OUT_DIR"

        apptainer exec \
            --bind "$DB_PATH":"$DB_PATH" \
            --bind "$OUT_DIR":"$OUT_DIR" \
            {params.container} \
            funannotate annotate \
                -i "$OUT_DIR" \
                --database "$DB_PATH" \
                --cpus {threads} \
                --busco_db {params.busco_db} \
                --eggnog "$EGGNOG" \
                --iprscan "$IPRSCAN" \
                --signalp "$SIGNALP" \
                --force

        out_path="$OUT_DIR/annotate_results"

        real_gff=$(find "$out_path" -maxdepth 1 -name "*.gff3" | head -n1 || true)
        real_gbk=$(find "$out_path" -maxdepth 1 -name "*.gbk" | head -n1 || true)
        real_prot=$(find "$out_path" -maxdepth 1 -name "*.proteins.fa" | head -n1 || true)

        if [ -n "$real_gff" ]; then
            ln -sfn "$(realpath "$real_gff")" {output.gff}
        else
            echo "ERROR: No GFF3 found in $out_path" >&2
            exit 1
        fi

        if [ -n "$real_gbk" ]; then
            ln -sfn "$(realpath "$real_gbk")" {output.gbk}
        else
            echo "ERROR: No GBK found in $out_path" >&2
            exit 1
        fi

        if [ -n "$real_prot" ]; then
            ln -sfn "$(realpath "$real_prot")" {output.proteins}
        else
            echo "ERROR: No proteins.fa found in $out_path" >&2
            exit 1
        fi

        cp -f "$DBCAN" "$out_path/{params.ref_name}.dbcan.overview.tsv"

        touch {output.done}
        """
        
# ============================================================================
# COMPARATIVE ANALYSIS
# ============================================================================
# busco_comparative: BUSCO completeness for the F. schrenkii and F. rosea assemblies
rule busco_comparative:
    input:
        assembly = lambda wc: COMPARATIVE[wc.species]["assembly"]
    output:
        summary = "results/comparative/busco/{species}/short_summary.txt",
        done = "results/comparative/busco/{species}/busco.done"
    params:
        lineage = config["qc"]["busco_lineage"],
        db_path = config["qc"]["busco_db_path"],
        out_path = "results/comparative/busco/{species}",
        out_name = "{species}"
    container:
        config["containers"]["busco"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.out_path}"
        ASSEMBLY=$(realpath {input.assembly})

        busco \
            -i "$ASSEMBLY" \
            -l "{params.lineage}" \
            -o "{params.out_name}" \
            --out_path "{params.out_path}" \
            --download_path "{params.db_path}" \
            -m genome \
            --cpu {threads} \
            --force

        mkdir -p "$(dirname {output.summary})"
        find "{params.out_path}/{params.out_name}" -name "short_summary*.txt" -exec cp {{}} "{output.summary}" \;

        test -s "{output.summary}"
        touch "{output.done}"
        """

# comparative_assembly_summary: combined assembly stats and BUSCO across all three species
rule comparative_assembly_summary:
    input:
        pin_assembly = COMPARATIVE["F_pinicola"]["assembly"],
        sch_assembly = COMPARATIVE["F_schrenkii"]["assembly"],
        ros_assembly = COMPARATIVE["F_rosea"]["assembly"],
        pin_busco = COMPARATIVE["F_pinicola"]["busco_summary"],
        sch_busco = COMPARATIVE["F_schrenkii"]["busco_summary"],
        ros_busco = COMPARATIVE["F_rosea"]["busco_summary"]
    output:
        table = "results/comparative/tables/comparative_assembly_summary.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/comparative_assembly_summary.py"
        
# species_specific_ipr_summary: InterPro domain summary for species-specific orthogroups
rule species_specific_ipr_summary:
    input:
        genes = "results/comparative/tables/species_specific_orthogroups.tsv",
        pin_ipr = f"results/annotation/{REFERENCE}/interproscan/interproscan.tsv",
        sch_ipr = "comparative/F_schrenkii/Fompi3_GeneCatalog_proteins_20120705_IPR.tab",
        ros_ipr = "comparative/F_rosea/Fomro1_GeneCatalog_proteins_20180323_IPR.tab"
    output:
        table = "results/comparative/tables/species_specific_ipr_summary.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/species_specific_ipr_summary.py"

# dbcan_comparative: dbCAN CAZyme annotation for F. schrenkii and F. rosea
rule dbcan_comparative:
    input:
        proteins = lambda wc: COMPARATIVE[wc.species]["proteins"]
    output:
        overview = "results/comparative/dbcan/{species}/overview.tsv",
        hmm = "results/comparative/dbcan/{species}/dbCAN_hmm_results.tsv",
        sub = "results/comparative/dbcan/{species}/dbCANsub_hmm_results.tsv",
        diamond = "results/comparative/dbcan/{species}/diamond.out",
        done = "results/comparative/dbcan/{species}/dbcan.done"
    params:
        outdir = "results/comparative/dbcan/{species}",
        db_dir = config["annotation"]["dbcan"]["db_dir"],
        evalue = config["annotation"]["dbcan"]["evalue"],
        coverage = config["annotation"]["dbcan"]["coverage"]
    threads: config["threads"]["medium"]
    log:
        "logs/dbcan/{species}.log"
    container:
        config["containers"]["dbcan"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        mkdir -p "$(dirname {log})"

        PROTEINS=$(realpath {input.proteins})
        OUTDIR=$(realpath {params.outdir})
        DB_DIR=$(realpath {params.db_dir})

        run_dbcan CAZyme_annotation \
            --input_raw_data "$PROTEINS" \
            --output_dir "$OUTDIR" \
            --db_dir "$DB_DIR" \
            --mode protein \
            --methods hmm --methods diamond --methods dbCANsub \
            --e_value_threshold_dbcan {params.evalue} \
            --coverage_threshold_dbcan {params.coverage} \
            --e_value_threshold_dbsub {params.evalue} \
            --coverage_threshold_dbsub {params.coverage} \
            --threads {threads} \
            2>&1 | tee {log}

        test -f {output.overview}
        touch {output.done}
        """

# dbcan_class_summary: per-class CAZyme counts across all three species
rule dbcan_class_summary:
    input:
        pinicola = f"results/annotation/{REFERENCE}/dbcan/overview.tsv",
        schrenkii = "results/comparative/dbcan/F_schrenkii/overview.tsv",
        rosea = "results/comparative/dbcan/F_rosea/overview.tsv"
    output:
        summary = "results/comparative/tables/dbcan_class_counts.tsv"
    params:
        species = {
            "F_pinicola": f"results/annotation/{REFERENCE}/dbcan/overview.tsv",
            "F_schrenkii": "results/comparative/dbcan/F_schrenkii/overview.tsv",
            "F_rosea": "results/comparative/dbcan/F_rosea/overview.tsv"
        }
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/dbcan_class_summary.py"

# dbcan_family_summary: per-family CAZyme counts across all three species
rule dbcan_family_summary:
    input:
        pinicola = f"results/annotation/{REFERENCE}/dbcan/overview.tsv",
        schrenkii = "results/comparative/dbcan/F_schrenkii/overview.tsv",
        rosea = "results/comparative/dbcan/F_rosea/overview.tsv"
    output:
        summary = "results/comparative/tables/dbcan_family_counts.tsv"
    params:
        species = {
            "F_pinicola": f"results/annotation/{REFERENCE}/dbcan/overview.tsv",
            "F_schrenkii": "results/comparative/dbcan/F_schrenkii/overview.tsv",
            "F_rosea": "results/comparative/dbcan/F_rosea/overview.tsv"
        }
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/dbcan_family_summary.py"

# plot_dbcan_family_heatmap: heatmap of the most variable CAZy families
rule plot_dbcan_family_heatmap:
    input:
        table = "results/comparative/tables/dbcan_family_counts.tsv"
    output:
        plot = "results/comparative/plots/dbcan_family_heatmap.png"
    conda:
        "envs/r_plotting.yaml"
    script:
        "scripts/plot_dbcan_family_heatmap.R"


# signalp_comparative: signal peptide prediction for F. schrenkii and F. rosea
rule signalp_comparative:
    input:
        proteins = lambda wc: COMPARATIVE[wc.species]["proteins"]
    output:
        summary = "results/comparative/signalp/{species}/{species}_signalp_summary.txt",
        done = "results/comparative/signalp/{species}/signalp.done"
    params:
        outdir = "results/comparative/signalp/{species}",
        signalp_dir = config["annotation"]["signalp"]["dir"],
        organism = "euk"
    threads: 4
    log:
        "logs/signalp_comparative/{species}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname {log})"

        OUTDIR=$(realpath "{params.outdir}")
        PROTEINS=$(realpath "{input.proteins}")
        SIGNALP_DIR=$(realpath "{params.signalp_dir}")
        LOGFILE=$(realpath "{log}")
        SUMMARY_OUT=$(realpath "$(dirname {output.summary})")/$(basename {output.summary})
        DONE_OUT=$(realpath "$(dirname {output.done})")/$(basename {output.done})

        cd "$SIGNALP_DIR/bin"

        ./signalp \
            -fasta "$PROTEINS" \
            -org {params.organism} \
            -format short \
            -prefix "$OUTDIR/{wildcards.species}_signalp" \
            2>&1 | tee "$LOGFILE"

        if [ -f "$OUTDIR/{wildcards.species}_signalp_summary.signalp5" ]; then
            cp "$OUTDIR/{wildcards.species}_signalp_summary.signalp5" "$SUMMARY_OUT"
        else
            echo "ERROR: Expected SignalP summary file not found in $OUTDIR" >&2
            ls -lah "$OUTDIR" >&2 || true
            exit 1
        fi

        test -s "$SUMMARY_OUT"
        touch "$DONE_OUT"
        """

# merge_secreted_cazymes: intersect dbCAN consensus with signal peptides per species
rule merge_secreted_cazymes:
    input:
        dbcan = dbcan_overview_for_species,
        signalp = signalp_summary_for_species
    output:
        secreted_cazymes = "results/comparative/secreted_cazymes_{species}.tsv"
    log:
        "logs/merge_secreted_cazymes/{species}.log"
    shell:
        r"""
        set -euo pipefail
        python scripts/merge_secreted_cazymes.py {input.dbcan} {input.signalp} {output.secreted_cazymes} > {log} 2>&1
        """

# secreted_dbcan_family_summary: per-family secreted CAZyme counts across species
rule secreted_dbcan_family_summary:
    input:
        pinicola = "results/comparative/secreted_cazymes_F_pinicola.tsv",
        schrenkii = "results/comparative/secreted_cazymes_F_schrenkii.tsv",
        rosea = "results/comparative/secreted_cazymes_F_rosea.tsv"
    output:
        summary = "results/comparative/tables/secreted_dbcan_family_counts.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/secreted_dbcan_summary.py"

# plot_secreted_dbcan_family_heatmap: heatmap of the most variable secreted CAZy families
rule plot_secreted_dbcan_family_heatmap:
    input:
        table = "results/comparative/tables/secreted_dbcan_family_counts.tsv"
    output:
        plot = "results/comparative/plots/secreted_dbcan_family_heatmap.png"
    conda:
        "envs/r_plotting.yaml"
    script:
        "scripts/plot_secreted_dbcan_heatmap.R"

# prep_orthofinder: stage species protein FASTA files for OrthoFinder
rule prep_orthofinder:
    input:
        pinicola = COMPARATIVE["F_pinicola"]["proteins"],
        schrenkii = COMPARATIVE["F_schrenkii"]["proteins"],
        rosea = COMPARATIVE["F_rosea"]["proteins"]
    output:
        outdir = directory("results/comparative/orthofinder_input")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.outdir}

        cp "$(realpath {input.pinicola})"  {output.outdir}/F_pinicola.fa
        cp "$(realpath {input.schrenkii})" {output.outdir}/F_schrenkii.fa
        cp "$(realpath {input.rosea})"     {output.outdir}/F_rosea.fa
        """

# orthofinder: orthology inference across all three species
rule orthofinder:
    input:
        indir = "results/comparative/orthofinder_input"
    output:
        results_dir = directory("results/comparative/orthofinder_results"),
        orthogroups = "results/comparative/orthofinder_results/Orthogroups/Orthogroups.tsv",
        gene_counts = "results/comparative/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv",
        single_copy = "results/comparative/orthofinder_results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
        stats = "results/comparative/orthofinder_results/Comparative_Genomics_Statistics/Statistics_Overall.tsv",
        done = "results/comparative/orthofinder.done"
    threads: config["threads"]["heavy"]
    container: config["containers"]["orthofinder"]
    shell:
        r"""
        set -euo pipefail

        orthofinder -f {input.indir} -t {threads} -a {threads}

        LATEST=$(ls -dt {input.indir}/OrthoFinder/Results_* | head -n1)
        [ -n "$LATEST" ] || (echo "ERROR: No OrthoFinder Results_* directory found" && exit 1)

        rm -rf results/comparative/orthofinder_results
        mkdir -p results/comparative
        cp -r "$LATEST" results/comparative/orthofinder_results

        test -f results/comparative/orthofinder_results/Orthogroups/Orthogroups.tsv
        test -f results/comparative/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv
        test -f results/comparative/orthofinder_results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
        test -f results/comparative/orthofinder_results/Comparative_Genomics_Statistics/Statistics_Overall.tsv

        touch {output.done}
        """

# orthofinder_presence_matrix: convert orthogroups to presence/absence matrix
rule orthofinder_presence_matrix:
    input:
        gene_counts = "results/comparative/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv"
    output:
        presence = "results/comparative/tables/orthogroups_presence.tsv"
    params:
        species = ["F_pinicola", "F_schrenkii", "F_rosea"]
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/orthofinder_presence_matrix.py"

# orthogroup_intersection_summary: count orthogroups in each species set intersection
rule orthogroup_intersection_summary:
    input:
        presence = "results/comparative/tables/orthogroups_presence.tsv"
    output:
        summary = "results/comparative/tables/orthogroup_intersection_summary.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/orthogroup_intersection_summary.py"

# plot_orthogroup_upset: UpSet plot of orthogroup sharing across species
rule plot_orthogroup_upset:
    input:
        presence = "results/comparative/tables/orthogroups_presence.tsv"
    output:
        plot = "results/comparative/plots/orthogroups_upset.png"
    conda:
        "envs/r_plotting.yaml"
    script:
        "scripts/plot_orthogroup_upset.R"

# orthogroup_species_specific_genes: extract genes from species-unique orthogroups
rule orthogroup_species_specific_genes:
    input:
        orthogroups = "results/comparative/orthofinder_results/Orthogroups/Orthogroups.tsv",
        presence = "results/comparative/tables/orthogroups_presence.tsv"
    output:
        table = "results/comparative/tables/species_specific_orthogroups.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/orthogroup_species_specific_genes.py"

# orthogroup_summary_table: combined orthogroup summary statistics
rule orthogroup_summary_table:
    input:
        presence = "results/comparative/tables/orthogroups_presence.tsv",
        single_copy = "results/comparative/orthofinder_results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
    output:
        table = "results/comparative/tables/orthogroup_summary_table.tsv"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/orthogroup_summary_table.py"

# eggnog_map_comparative: eggNOG-mapper for F. schrenkii and F. rosea
rule eggnog_map_comparative:
    input:
        faa = lambda wc: COMPARATIVE[wc.species]["proteins"]
    output:
        ann = "results/comparative/eggnog/{species}/annotations.emapper.annotations"
    params:
        data_dir = config["annotation"]["eggnog"]["db_dir"],
        outdir = "results/comparative/eggnog/{species}",
        prefix = "annotations",
        tax = config["annotation"]["eggnog"]["tax_scope"]
    threads: config["threads"]["medium"]
    container: config["containers"]["eggnog-mapper"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        emapper.py \
            -i {input.faa} \
            --itype proteins \
            -o {params.prefix} \
            --output_dir {params.outdir} \
            --cpu {threads} \
            --usemem \
            --tax_scope {params.tax} \
            --data_dir "{params.data_dir}"
        """

# eggnog_cog_comparison: COG category counts across all three species
rule eggnog_cog_comparison:
    input:
        pin = f"results/annotation/{REFERENCE}/eggnog/annotations.emapper.annotations",
        sch = "results/comparative/eggnog/F_schrenkii/annotations.emapper.annotations",
        ros = "results/comparative/eggnog/F_rosea/annotations.emapper.annotations"
    output:
        table = "results/comparative/tables/cog_category_comparison.tsv",
        summary = "results/comparative/tables/cog_category_summary.txt"
    conda:
        "envs/plotting_python.yaml"
    script:
        "scripts/eggnog_cog_comparison.py"

# plot_cog_category_barplot: grouped barplot of COG category counts
rule plot_cog_category_barplot:
    input:
        table = "results/comparative/tables/cog_category_comparison.tsv"
    output:
        plot = "results/comparative/plots/cog_category_barplot.png"
    conda:
        "envs/r_plotting.yaml"
    script:
        "scripts/plot_cog_category_barplot.R"

# ============================================================================
# CIRCOS PLOT INPUT DATA GENERATION
# ============================================================================
# gene_density: per-window gene density from the funannotate GFF file
rule gene_density:
    input:
        gff3    = f"results/annotation/{REFERENCE}/annotate_results/{REFERENCE}.gff3",
        windows = config["external"]["genome_windows_10kb"]
    output:
        density = f"results/circos/{REFERENCE}/gene_density_10kb.tsv"
    container: config["containers"]["bedtools"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.density})

        # Extract gene features as BED
        awk -F'\t' '$3 == "gene" {{OFS="\t"; print $1, $4-1, $5}}' {input.gff3} \
            > {output.density}.tmp_genes.bed

        # Count genes per window
        bedtools coverage -a {input.windows} -b {output.density}.tmp_genes.bed \
            -counts \
            | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, $4}}' \
            > {output.density}

        rm -f {output.density}.tmp_genes.bed
        test -s {output.density}
        """

# parse_repeatmasker: convert RepeatMasker .out to BED with class labels
rule parse_repeatmasker:
    input:
        rm_out = f"results/annotation/{REFERENCE}/assembly.fasta.out"
    output:
        bed = f"results/circos/{REFERENCE}/repeats_by_class.bed"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed})

        python3 scripts/parse_repeatmasker_out.py \
            --input {input.rm_out} \
            --output {output.bed}

        test -s {output.bed}
        """

# repeat_bp_density: per-window bp coverage by repeat class for the circos plot
rule repeat_bp_density:
    input:
        repeats = f"results/circos/{REFERENCE}/repeats_by_class.bed",
        windows = config["external"]["genome_windows_10kb"]
    output:
        density = f"results/circos/{REFERENCE}/repeat_bp_density_10kb.tsv"
    container: config["containers"]["bedtools"]
    shell:
        r"""
        set -euo pipefail

        echo -e "chrom\tstart\tend\tLTR\tDNA_TE\tLINE\tUnclassified\tOther" \
            > {output.density}

        # Split repeats by class
        for CLASS in LTR DNA_TE LINE Unclassified Other; do
            awk -v c="$CLASS" '$4 == c' {input.repeats} \
                > {output.density}.tmp_${{CLASS}}.bed || true
        done

        # For the first class keep genomic coordinates + covered bp
        bedtools coverage -a {input.windows} \
            -b {output.density}.tmp_LTR.bed \
            | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, $5}}' \
            > {output.density}.tmp_cov_LTR

        # For remaining classes keep only covered bp column
        for CLASS in DNA_TE LINE Unclassified Other; do
            bedtools coverage -a {input.windows} \
                -b {output.density}.tmp_${{CLASS}}.bed \
                | awk 'BEGIN{{OFS="\t"}} {{print $5}}' \
                > {output.density}.tmp_cov_${{CLASS}}
        done

        paste {output.density}.tmp_cov_LTR \
              {output.density}.tmp_cov_DNA_TE \
              {output.density}.tmp_cov_LINE \
              {output.density}.tmp_cov_Unclassified \
              {output.density}.tmp_cov_Other \
            >> {output.density}

        rm -f {output.density}.tmp_*
        test -s {output.density}
        """
