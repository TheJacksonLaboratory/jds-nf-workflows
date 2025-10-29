import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
HAPLOTYPE RECONSTRUCTION PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--rerun                         ${params.rerun}
--correct_ids                   ${params.correct_ids}
--remove_markers                ${params.remove_markers}
"""
}