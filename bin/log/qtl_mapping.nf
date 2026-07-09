def param_log(){

def message = ""

message =  """
QTL MAPPING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--n_perms                       ${params.n_perms}
"""

log.info(message)
return(message)

}
