#!/usr/bin/env groovy

String javaBin = 'java'

def checkFailed = { file->
    def proc = "tail -n 2 ${file}".execute()
    proc.waitFor()
    return !proc.text.contains('End time:')
}//

if( !this.args || this.args.length!=1 ){
    println 'Get failed eqtl jobs listed in command file'
    println 'Usage: get_failed_eqtl_jobs.groovy <command file>'
    return 1
}

String commFile = this.args[0]
def writer = System.out

// extra option to add
def extraOption = '--isoforms=/home/users/ipb/lindo/filtered_isoforms/filtered_isoforms.final.txt'
def replace = ['-Xmx62g':'-Xmx20g']

new File(commFile).eachLine{ line->
    if( !line.startsWith(javaBin) ){
        writer.println(line)
    }
    else{
        def toks = line.split("\\s")
        def jobOutFile = new File(toks[toks.length-1])
        
        if( !jobOutFile.exists() || checkFailed(jobOutFile.absolutePath) ){
            if( extraOption ){
                int pos = line.indexOf('--props=')
                line = line.substring(0,pos)+extraOption+line.substring(pos-1)
            }
            
            if( replace ) {
                replace.each{k,v-> line = line.replace(k,v)}
            }
            
            writer.println(line)
        }
    }
}

writer.close()

return 0