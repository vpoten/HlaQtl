#!/usr/bin/env groovy

String CORREL_FILE = 'correlation.txt.filter'


if( !this.args || this.args.length!=3 ){
    println 'Search eQTLs in results folder for a SNPs or Isoforms list'
    println 'Usage: search_eqtls.groovy <dir> [--snps | --isoforms] <list>'
    return 1
}

String dir = this.args[0]
String mode = this.args[1]
String listFile = this.args[2]

def field = (mode=='--snps') ? 0 : 1

// fill snps/isoforms list
def list = [] as TreeSet
new File(listFile).eachLine{ list << it }

println "Start time: ${new Date()}\n"

new File(dir).eachDir{ fDir->
    def file = new File(fDir.absolutePath+'/'+CORREL_FILE)
    
    if( file.exists() ){
        def reader = file.newReader()
        reader.readLine()
        
        reader.splitEachLine("\\s"){ toks->
            if( toks[field] in list ){
                print file.absolutePath // print result dir name
                toks.each{ print '\t'; print it } // print correlation line
                print '\n'
            }
        }
        
        reader.close()
    }
}

println "End time: ${new Date()}"

return 0