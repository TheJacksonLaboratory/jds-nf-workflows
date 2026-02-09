import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
QTL MAPPING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--n_perms                       ${params.n_perms}
"""
}