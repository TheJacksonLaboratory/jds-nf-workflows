#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {AMPLICON_FINGERPRINT} from './workflows/amplicon_fingerprint'
include {AMPLICON_GENERIC} from './workflows/amplicon_generic'
include {ANCESTRY_RUN} from './workflows/ancestry'
include {ATAC} from './workflows/atac'
include {BAM_TO_FASTQ} from './supportworkflows/bam_to_fastq'
include {CHIPSEQ} from './workflows/chipseq'
include {CNV_ARRAY} from './workflows/cnv_array'
include {EMASE} from './workflows/emase'
include {GBRS} from './workflows/gbrs'
include {GENERATE_PSEUDOREFERENCE} from './workflows/generate_pseudoreference'
include {GENERATE_RNASEQ_INDEX} from './supportworkflows/generate_rnaseq_index'
include {GENERATE_RNASEQ_SIMREADS} from './supportworkflows/generate_rnaseq_simreads'
include {GENERATE_WES_SIMREADS} from './supportworkflows/generate_wes_simreads'
include {GENERATE_WGS_SIMREADS} from './supportworkflows/generate_wgs_simreads'
include {GERMLINE_SV} from "./workflows/germline_sv" // MMRSVDB
include {HAPLOTYPE_RECONSTRUCTION} from './workflows/haplotype_reconstruction'
include {JOINT_GVCF_CALLING} from './workflows/joint_gvcf_calling'
include {LCWGS_HR} from './workflows/lcwgs_hr'
include {MITOCHONDRIA_VARIANT_CALLING} from './workflows/mitochondria_variant_calling'
include {PREPARE_EMASE} from './workflows/prepare_emase'
include {PREP_DO_GBRS_INPUT} from './supportworkflows/prep_do_gbrs_inputs'
include {PTA} from './workflows/pta'
include {QTL_MAPPING} from './workflows/qtl_mapping'
include {REANNOTATE_PTA} from './supportworkflows/reannotate_pta'
include {RNA_FUSION} from './workflows/rna_fusion'
include {RNASEQ} from './workflows/rnaseq'
include {RRBS} from './workflows/rrbs'
include {SMRNASEQ} from './workflows/smrnaseq'
include {SOMATIC_WES} from './workflows/somatic_wes'
include {SOMATIC_WES_PTA} from './workflows/somatic_wes_pta'
include {WES} from './workflows/wes'
include {WGS} from './workflows/wgs'
include {WGS_LONG_READ} from './workflows/wgs_long_read'
include {WGS_SV_BAM} from './workflows/wgs_sv_bam'

// conditional to launch appropriate workflow
workflow{
  println([
    '',
    '_____ _____  ______     __   _ _____     _  _  _  ____   _____ _     _ _____        ____  _  _  _ ______',
    '  |   |    \\ |_____ ___ | \\  | |____ ___ |  |  | |    | |____/ |____/  |____ |     |    | |  |  | |_____',
    '__|   |____/ _____|     |  \\_| |         |__|__| |____| |   \\_ |    \\_ |     |____ |____| |__|__| _____|'
  ].join('\n'))

  // The logo looks incorrect above, but the backslashes need to be escaped in the string. The actual string is a logo for the pipeline.

  if (params.workflow == "amplicon" || params.workflow == "amplicon_fingerprint"){
     AMPLICON_FINGERPRINT()
  }
  if (params.workflow == "amplicon_generic"){
    AMPLICON_GENERIC()
  }
  if (params.workflow == "ancestry"){
    ANCESTRY_RUN()
  }
  if (params.workflow == "atac"){
    ATAC()
  }
  if (params.workflow == "bam_to_fastq"){
    BAM_TO_FASTQ()
  }
  if (params.workflow == "chipseq"){
    CHIPSEQ()
  }
  if (params.workflow == "cnv_array"){
    CNV_ARRAY()
  }
  if (params.workflow == "emase"){
    EMASE()
  }
  if (params.workflow == "gbrs"){
    GBRS()
  }
  if (params.workflow == "generate_pseudoreference") {
    GENERATE_PSEUDOREFERENCE()
  }
  if (params.workflow == "generate_rnaseq_index"){
    GENERATE_RNASEQ_INDEX()
  }
  if (params.workflow == "generate_rnaseq_simreads"){
    GENERATE_RNASEQ_SIMREADS()
  }
  if (params.workflow == "generate_wes_simreads"){
    GENERATE_WES_SIMREADS()
  }
  if (params.workflow == "generate_wgs_simreads"){
    GENERATE_WGS_SIMREADS()
  }
  if (params.workflow == "germline_sv"){
    GERMLINE_SV()
  }
  if (params.workflow == "haplotype_reconstruction"){
    HAPLOTYPE_RECONSTRUCTION()
  }
  if (params.workflow == "joint_gvcf_calling"){
    JOINT_GVCF_CALLING()
  }
  if (params.workflow == "lcwgs_hr"){
    LCWGS_HR()
  }
  if (params.workflow == "mitochondria_variant_calling"){
    MITOCHONDRIA_VARIANT_CALLING()
  }
  if (params.workflow == "prepare_emase"){
    PREPARE_EMASE()
  }
  if (params.workflow == "prep_do_gbrs_inputs"){
    PREP_DO_GBRS_INPUT()
  }
  if (params.workflow == "pta"){
    PTA()
  }
  if (params.workflow == "qtl_mapping"){
    QTL_MAPPING()
  }
  if (params.workflow == "reannotate_pta"){
    REANNOTATE_PTA()
  }
  if (params.workflow == "rna_fusion"){
    RNA_FUSION()
  }
  if (params.workflow == "rnaseq"){
    RNASEQ()
  }
  if (params.workflow == "rrbs"){
    RRBS()
  }
  if (params.workflow == "smrnaseq"){
    SMRNASEQ()
  }
  if (params.workflow == "wes"){
    WES()
  }
  if (params.workflow == "somatic_wes"){
    SOMATIC_WES()
  }
  if (params.workflow == "somatic_wes_pta"){
    SOMATIC_WES_PTA()
  }
  if (params.workflow == "wgs"){
    WGS()
  }
  if (params.workflow == "wgs_long_read"){
    WGS_LONG_READ()
  }
  if (params.workflow == "wgs_sv_bam"){
    WGS_SV_BAM()
  }
  else {
    exit 1, "Invalid workflow specified: ${params.workflow}. Please specify a valid workflow using the `--workflow` parameter. See `--help` for information."
  }
}
