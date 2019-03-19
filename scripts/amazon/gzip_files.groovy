#!/usr/bin/env groovy

if( !this.args || this.args.length!=1 ){
    println 'Gzip matching files in directory'
    println 'usage: gzip_files.groovy <directory>'
    return 1
}

def patterns = [/correlation\.txt/, /debug_corrData\.txt/, /iso_merged\.tracking/]

String dir = this.args[0]

    
if( !dir.endsWith('/') )
    dir += '/'
    

new File(dir).eachFileRecurse{ file->
    if( file.file && patterns.any{file.name==~it} ){
        println "gzip ${file.absolutePath}"
        "gzip -f ${file.absolutePath}".execute().waitFor()
    }
}

return 0