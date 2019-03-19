#!/usr/bin/env groovy

String CORREL_FILE = 'correlation.txt.filter'

def buildKey = { a, b->
    "${a}_${b}"
}//

if( !this.args || this.args.length!=2 ){
    println 'Join eQTLs in results folder'
    println 'Usage: join_eqtls.groovy <dir> <out_file>'
    return 1
}

String dir = this.args[0]
String outFile = this.args[1]


println "Start time: ${new Date()}\n"

def writer = new PrintWriter(outFile)
def written = [] as TreeSet
int count = 0

new File(dir).eachDir{ fDir->
    def file = new File(fDir.absolutePath+'/'+CORREL_FILE)
    
    if( file.exists() ){
        count++
        def reader = file.newReader()
        reader.readLine() //skip header
        
        reader.eachLine{ line->
            def toks = line.split("\\s", 3)
            def key = buildKey(toks[0], toks[1])
            
            if( !(key in written) ){
                writer.println(line)
                written << key
            }
        }
        
        reader.close()
    }
}

writer.close()

println "${count} files processed."
println "${written.size()} eQTLs written."

println "End time: ${new Date()}"

return 0