import Logos

logo = new Logo()
println '\n'
println logo.show()

if (params.gen_org != "human" && params.gen_org != "mouse") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', or 'human'" 
}

def param_log(){

if (params.csv_input)
log.info """
EMASE RUN PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--sample_folder                 ${params.sample_folder}

--csv_input                     ${params.csv_input}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else
log.info """
EMASE RUN PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}


Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}