#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Microscape Amplicon Pipeline - Nextflow DSL2
// ============================================================================
//
// Amplicon sequencing analysis pipeline: demultiplexing, primer removal,
// DADA2 denoising, chimera removal, taxonomy assignment, phylogeny,
// normalization, clustering, network analysis, and Shiny visualization.
//
// Usage:
//   nextflow run main.nf --input /path/to/reads --primers primers-all.fa -resume
//
// ============================================================================

// ============================================================================
// Help message
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Microscape Amplicon Pipeline
     https://github.com/rec3141/microscape
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads [options] -resume

    Required:
      --input DIR          Input directory containing demultiplexed paired-end
                           *.fastq.gz files (PLATE_WELL_R{1,2}.fastq.gz)

    Primer removal:
      --primer_auto        Auto-select primer file by filename prefix [default: true]
      --primers PATH       Primer FASTA (overrides auto-selection)
      --primer_error_rate  Cutadapt error rate [default: 0.12]

    Demultiplexing (optional):
      --run_demultiplex    Enable demultiplexing with Mr_Demuxy [default: false]
      --forward_bcs PATH   Forward barcode file
      --reverse_bcs PATH   Reverse barcode file

    DADA2 parameters:
      --maxEE N            Max expected errors [default: 2]
      --truncQ N           Truncation quality [default: 11]
      --maxN N             Max Ns [default: 0]
      --truncLen_fwd N     Truncate forward reads at position N [default: 0 = off]
      --truncLen_rev N     Truncate reverse reads at position N [default: 0 = off]
      --min_overlap N      Min overlap for pair merging [default: 10]

    QC filtering:
      --min_seq_length N   Min sequence length [default: 50]
      --min_reads N        Min reads per sample [default: 100]
      --min_samples N      Min samples per ESV [default: 2]
      --min_seqs N         Min total reads per ESV [default: 3]

    Output:
      --outdir DIR         Output directory [default: results]
      --store_dir DIR      Persistent cache directory [default: off]

    Taxonomy:
      --ref_databases STR  Reference databases (semicolon-separated):
                           "name:path:Level1,Level2,..." [default: none]
      --run_phylogeny      Build phylogenetic tree [default: false]

    Metadata:
      --metadata PATH      MIMARKS-compliant TSV/CSV metadata file [default: none]
      --sample_id_column S Column name matching sample IDs [default: sample_name]

    Network:
      --min_prevalence N   Min samples per ASV for SparCC [default: 5]

    Language:
      --lang STR           'R' or 'python' [default: R]
      --dada2_path PATH Path to dada2 Python package (rec3141/dada2 gpu-python branch) [default: none]

    Resources:
      --threads N          General thread count [default: 8]
      --dada_cpus N        CPUs for DADA2 processes [default: 8]
      --dada_memory S      Memory for DADA2 processes [default: '16 GB']
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

// ============================================================================
// Parameter validation
// ============================================================================

if (!params.input) {
    log.error "ERROR: --input is required. Provide path to directory containing *.fastq.gz files. Run with --help for usage."
    System.exit(1)
}

if (params.run_demultiplex && (!params.forward_bcs || !params.reverse_bcs)) {
    log.error "ERROR: --forward_bcs and --reverse_bcs are required when --run_demultiplex is enabled."
    System.exit(1)
}

// ============================================================================
// Module imports
// ============================================================================

// Stage A: Preprocessing and DADA2
include { DEMULTIPLEX }       from './modules/demultiplex'
include { REMOVE_PRIMERS }    from './modules/primers'
include { DADA2_FILTER_TRIM } from './modules/dada2'
include { DADA2_LEARN_ERRORS } from './modules/dada2'
include { DADA2_DENOISE }     from './modules/dada2'
include { MERGE_SEQTABS }     from './modules/merge'
include { REMOVE_CHIMERAS }   from './modules/merge'
include { FILTER_SEQTAB }     from './modules/merge'

// Stage B: Taxonomy, phylogeny, renormalization
include { ASSIGN_TAXONOMY }   from './modules/taxonomy'
include { BUILD_PHYLOGENY }   from './modules/phylogeny'
include { RENORMALIZE }       from './modules/renormalize'

// Stage C: Metadata, clustering, network, visualization
include { LOAD_METADATA }     from './modules/metadata'
include { CLUSTER_TSNE }      from './modules/cluster'
include { NETWORK_SPARCC }    from './modules/network'
include { BUILD_VIZ }         from './modules/shiny'

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    // 1. Discover input reads
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
    }

    // 2. Optional demultiplexing
    if (params.run_demultiplex) {
        ch_raw = Channel.fromPath("${params.input}/*.fastq.gz")
            .map { fastq ->
                def name = fastq.baseName.replace('.fastq', '')
                def lane = name.split('_')[0]
                [[id: lane], fastq]
            }
            .groupTuple()

        DEMULTIPLEX(ch_raw,
                    file(params.forward_bcs),
                    file(params.reverse_bcs))
        ch_demuxed = DEMULTIPLEX.out.reads.flatten()
            .filter { it.name.endsWith('.fastq.gz') }
            .map { fastq ->
                def name = fastq.baseName.replace('.fastq', '').replace('.ctrimmed', '')
                [[id: name], fastq]
            }
    } else {
        // Pair up R1/R2 files by sample prefix
        ch_demuxed = Channel.fromFilePairs("${params.input}/*_R{1,2}*.fastq.gz", flat: true)
            .map { sample_id, r1, r2 ->
                [[id: sample_id], r1, r2]
            }
    }

    // 3. Remove primers with cutadapt (skip if input is already trimmed)
    if (params.skip_primer_removal) {
        ch_trimmed = ch_demuxed
            .map { meta, r1, r2 ->
                def plate = meta.id.split('_')[0]
                def new_meta = meta + [plate: plate]
                [new_meta, r1, r2]
            }
    } else {
        REMOVE_PRIMERS(ch_demuxed)
        ch_trimmed = REMOVE_PRIMERS.out.reads
            .map { meta, r1, r2 ->
                def plate = meta.id.split('_')[0]
                def new_meta = meta + [plate: plate]
                [new_meta, r1, r2]
            }
    }

    // 4a. Filter and trim per-sample
    DADA2_FILTER_TRIM(ch_trimmed)

    // 4b. Group filtered reads by plate
    ch_by_plate = DADA2_FILTER_TRIM.out.reads
        .filter { meta, r1, r2 -> r1.size() > 0 && r2.size() > 0 }
        .map { meta, r1, r2 -> [meta.plate, r1, r2] }
        .groupTuple(by: 0)
        .map { plate, r1s, r2s ->
            [[id: plate, plate: plate], r1s, r2s]
        }

    // 4c. Learn error rates per-plate
    DADA2_LEARN_ERRORS(ch_by_plate)

    // 4d. Denoise per-plate using learned error rates
    //     Join filtered reads with error models by plate name
    ch_plate_reads = ch_by_plate
        .map { meta, r1s, r2s -> [meta.plate, meta, r1s, r2s] }

    ch_plate_errors = DADA2_LEARN_ERRORS.out.error_models_rds
        .map { plate, meta, errF, errR -> [plate, errF, errR] }

    ch_denoise_input = ch_plate_reads
        .join(ch_plate_errors, by: 0)
        .map { plate, meta, r1s, r2s, errF, errR ->
            [meta, r1s, r2s, errF, errR]
        }

    DADA2_DENOISE(ch_denoise_input)

    // 5. Merge all sequence tables and remove chimeras
    ch_all_seqtabs = DADA2_DENOISE.out.seqtab.map { meta, rds -> rds }.collect()

    MERGE_SEQTABS(ch_all_seqtabs)

    // Per-plate chimera removal already happens in DADA2_DENOISE (both R and Python).
    // Post-merge chimera pass catches cross-plate chimeras.
    def effectiveEngine = params.dada_engine ?: params.lang
    if (effectiveEngine == 'python') {
        REMOVE_CHIMERAS(MERGE_SEQTABS.out.seqtab)
        FILTER_SEQTAB(REMOVE_CHIMERAS.out.seqtab)
    } else {
        FILTER_SEQTAB(MERGE_SEQTABS.out.seqtab)
    }

    // ======================================================================
    // Stage B: Taxonomy, phylogeny, renormalization
    // ======================================================================

    // 6. Assign taxonomy — one task per reference database (parallel)
    //    Build a channel of (db_name, db_path, tax_levels) tuples from params
    if (params.ref_databases) {
        def db_entries = params.ref_databases.tokenize(';')
        ch_databases = Channel.from(db_entries).map { entry ->
            // Format: "name:path:Level1,Level2,..." or "name:path" (no levels)
            def parts = entry.tokenize(':')
            def db_name = parts[0].trim()
            def db_path = file(parts[1].trim())
            def tax_levels = parts.size() > 2 ? parts[2].trim() : null
            [db_name, db_path, tax_levels]
        }

        ASSIGN_TAXONOMY(ch_databases, FILTER_SEQTAB.out.seqtab)

        // 9. Renormalize using the primary taxonomy database
        //    Use the first database as the primary for group classification
        ch_primary_tax = ASSIGN_TAXONOMY.out.taxonomy.first()
        RENORMALIZE(FILTER_SEQTAB.out.seqtab, ch_primary_tax)
    }

    // 7. Build phylogeny (optional)
    if (params.run_phylogeny) {
        BUILD_PHYLOGENY(FILTER_SEQTAB.out.seqtab)
    }

    // ======================================================================
    // Stage C: Metadata, clustering, network, visualization
    // ======================================================================

    // 8. Load metadata (optional — if metadata file provided)
    if (params.metadata) {
        LOAD_METADATA(FILTER_SEQTAB.out.seqtab,
                      file(params.metadata),
                      params.sample_id_column)
        ch_metadata = LOAD_METADATA.out.metadata
    } else {
        ch_metadata = Channel.empty()
    }

    // 10. t-SNE clustering (samples + ASVs)
    CLUSTER_TSNE(FILTER_SEQTAB.out.seqtab)

    // 11. SparCC network analysis (requires renormalized data)
    if (params.ref_databases) {
        NETWORK_SPARCC(RENORMALIZE.out.merged, params.min_prevalence)

        // 12. Build Shiny app (bundle all outputs)
        //     Collect all taxonomy files into a directory
        ch_tax_files = ASSIGN_TAXONOMY.out.taxonomy
            .map { db_name, tax, boot -> [tax, boot] }
            .flatten()
            .collect()

        BUILD_VIZ(
            FILTER_SEQTAB.out.seqtab,
            RENORMALIZE.out.merged,
            ch_tax_files,
            ch_metadata.ifEmpty(file('NO_METADATA')),
            CLUSTER_TSNE.out.sample_tsne,
            CLUSTER_TSNE.out.seq_tsne,
            NETWORK_SPARCC.out.correlations
        )
    }
}

// ============================================================================
// Pipeline completion handler
// ============================================================================

workflow.onComplete {
    def msg = """\
        Pipeline completed at : ${workflow.complete}
        Duration              : ${workflow.duration}
        Success               : ${workflow.success}
        Exit status           : ${workflow.exitStatus}
        Output directory      : ${params.outdir}
        """.stripIndent()

    println msg

    if (!workflow.success) {
        println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
    }
}

workflow.onError {
    println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
}
