#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}
else if (params.workflow == "wes"){
  include {WES} from './workflows/wes'
}
else if (params.workflow == "somatic_wes"){
  include {SOMATIC_WES} from './workflows/somatic_wes'
}
else if (params.workflow == "somatic_wes_pta"){
  include {SOMATIC_WES_PTA} from './workflows/somatic_wes_pta'
}
else if (params.workflow == "wgs"){
  include {WGS} from './workflows/wgs'
}
else if (params.workflow == "rrbs"){
  include {RRBS} from './workflows/rrbs'
}
else if (params.workflow == "atac"){
  include {ATAC} from './workflows/atac'
}
else if (params.workflow == "chipseq"){
  include {CHIPSEQ} from './workflows/chipseq'
}
else if (params.workflow == "pta"){
  include {PTA} from './workflows/pta'
} 
else if (params.workflow == "rna_fusion"){
  include {RNA_FUSION} from './workflows/rna_fusion'
}
else if (params.workflow == "generate_pseudoreference"){
  include {GENERATE_PSEUDOREFERENCE} from './workflows/generate_pseudoreference'
}
else if (params.workflow == "prepare_emase"){
  include {PREPARE_EMASE} from './workflows/prepare_emase'
}
else if (params.workflow == "prep_do_gbrs_inputs"){
  include {PREP_DO_GBRS_INPUT} from './subworkflows/prep_do_gbrs_inputs'
}
else if (params.workflow == "emase"){
  include {EMASE} from './workflows/emase'
}
else if (params.workflow == "gbrs"){
  include {GBRS} from './workflows/gbrs'
}
else if (params.workflow == "amplicon"){
  include {AMPLICON} from './workflows/amplicon_fingerprint'
}
else if (params.workflow == "amplicon_generic"){
  include {AMPLICON} from './workflows/amplicon_generic'
}
else if (params.workflow == "ancestry"){
  include {ANCESTRY_RUN} from './workflows/ancestry'
}
else if (params.workflow == "germline_sv") { // MMRSVDB
  include {GERMLINE_SV} from "./workflows/germline_sv"
}
else if (params.workflow == "smrnaseq"){
  include {SMRNASEQ} from './workflows/smrnaseq'
}
else if (params.workflow == "cnv_array"){
  include {CNV_ARRAY} from './workflows/cnv_array'
}
else if (params.workflow == "generate_rnaseq_index"){
  include {GENERATE_RNASEQ_INDEX} from './subworkflows/generate_rnaseq_index'
}
else if (params.workflow == "joint_gvcf_calling"){
  include {JOINT_GVCF_CALLING} from './workflows/joint_gvcf_calling'
}
else if (params.workflow == "wgs_long_read"){
  include {wgs_long_read} from './workflows/wgs_long_read'
}
else if (params.workflow == "mitochondria_variant_calling"){
  include {MITOCHONDRIA_VARIANT_CALLING} from './workflows/mitochondria_variant_calling'
}
else if (params.workflow == "reannotate_pta"){
  include {REANNOTATE_PTA} from './subworkflows/reannotate_pta'
}
else if (params.workflow == "wgs_sv_bam"){
  include {WGS_SV_BAM} from './workflows/wgs_sv_bam'
}
else {
  // if workflow name is not supported: 
  exit 1, "ERROR: No valid pipeline called. '--workflow ${params.workflow}' is not a valid workflow name."
}

// conditional to launch appropriate workflow
workflow{
  if (params.workflow == "rnaseq"){
    RNASEQ()
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
  if (params.workflow == "rrbs"){
    RRBS()
    }
  if (params.workflow == "atac"){
    ATAC()
    } 
  if (params.workflow == "chipseq"){
    CHIPSEQ()
    }
  if (params.workflow == "pta"){
    PTA()
  } 
  if (params.workflow == "rna_fusion"){
    RNA_FUSION()
  }
  if (params.workflow == "generate_pseudoreference") {
    GENERATE_PSEUDOREFERENCE()
  }
  if (params.workflow == "prepare_emase"){
    PREPARE_EMASE()
  }
  if (params.workflow == "emase"){
    EMASE()
  }
  if (params.workflow == "gbrs"){
    GBRS()
  }
  if (params.workflow == "prep_do_gbrs_inputs"){
    PREP_DO_GBRS_INPUT()
  }
  if (params.workflow == "amplicon" || params.workflow == "amplicon_generic"){
    AMPLICON()
  }
  if (params.workflow == "ancestry"){
    ANCESTRY_RUN()
  }
  if (params.workflow == "germline_sv"){
    GERMLINE_SV()
  }
  if (params.workflow == "smrnaseq"){
    SMRNASEQ()
  }
  if (params.workflow == "cnv_array"){
    CNV_ARRAY()
  }
  if (params.workflow == "generate_rnaseq_index"){
    GENERATE_RNASEQ_INDEX()
  }
  if (params.workflow == "joint_gvcf_calling"){
    JOINT_GVCF_CALLING()
  }
  if (params.workflow == "wgs_long_read"){
    wgs_long_read()
  }
  if (params.workflow == "mitochondria_variant_calling"){
    MITOCHONDRIA_VARIANT_CALLING()
  }
  if (params.workflow == "reannotate_pta"){
    REANNOTATE_PTA()
  }
  if (params.workflow == "wgs_sv_bam"){
    WGS_SV_BAM()
  }
}
