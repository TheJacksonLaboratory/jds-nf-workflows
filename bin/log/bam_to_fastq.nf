def param_log(){

if (params.read_type != "PE" && params.read_type != "SE") {
    error "'--read_type': \"${params.read_type}\" is not valid, supported options are 'PE' or 'SE'" 
}


def message = ""

if (params.csv_input) {
message = """
GBRS RUN PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--read_type                     ${params.read_type}
--csv_input                     ${params.csv_input}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
} 

log.info(message)
return(message)

}

