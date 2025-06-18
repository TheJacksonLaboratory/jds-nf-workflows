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
--gen_org                       ${params.gen_org}
--n_perms                       ${params.n_perms}
"""
}