def final_run_report(message){
def report = """
==== Workflow Completion Report ====

Run Details:
===========================================
Run Name:         ${workflow.runName}
Session ID:       ${workflow.sessionId}

Command Line:     ${workflow.commandLine}

Script Name:      ${workflow.scriptName}
Script File:      ${workflow.scriptFile}
Script ID:        ${workflow.scriptId}
NextflowVer:      ${nextflow.version}

Output Dir:       ${workflow.outputDir}
Work Dir:         ${workflow.workDir}
Project Dir:      ${workflow.projectDir}
Profile:          ${workflow.profile}
Resume:           ${workflow.resume}

Run Time and Status Details:
===========================================
Start Time:       ${workflow.start}
Completion Time:  ${workflow.complete}
Run Duration:     ${workflow.duration}

Workflow Status:  ${workflow.success ? 'PASS' : 'FAIL'}
Exit Status:      ${workflow.exitStatus}
Error Message:    ${workflow.errorMessage ?: 'None'}
Error Report:     ${workflow.errorReport ?: 'None'}

GitHub Version Info (if available):
===========================================
Commit ID:        ${workflow.commitId ?: 'N/A'}
Revision:         ${workflow.revision ?: 'N/A'}
Repository:       ${workflow.repository ?: 'N/A'}

Manifest Info:
===========================================
Name:             ${workflow.manifest.name}
Repository:       ${workflow.manifest.homePage}
Description:      ${workflow.manifest.description}
Version:          ${workflow.manifest.version}
Contributors:     ${workflow.manifest.contributors.name.join(', ')}

Runtime Logged Info (shown to stdout during run):
===========================================
${message}

"""
def fileName = "workflow_report_${workflow.runName}_${workflow.sessionId}.txt"
new File(fileName).text = report
}
