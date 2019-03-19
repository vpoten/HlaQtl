#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

String CORREL_FILE = 'correlation.txt.gz'

def hash = { string->
    long h = 1125899906842597L; // prime
    int len = string.length();

    for (int i = 0; i < len; i++) {
        h = 31*h + string.charAt(i);
    }
    return h;
}//
    
def buildKey = { a, b->
    hash("${a}_${b}")
}//

if( !this.args || this.args.length!=4 ){
    println 'Join all eQTLs (significant and non-significant) in results folder'
    println 'Usage: join_all_eqtls.groovy <eqtl dir> <isoform filter file> <type_label> <out_file>'
    return 1
}

String dir = this.args[0]
String isoFile = this.args[1]
String outFile = this.args[3]
String type = this.args[2]


println "Start time: ${new Date()}\n"

//get list of isoforms
def isoSet = [] as TreeSet
def reader = new File(isoFile).newReader()

reader.readLine()//skip header
reader.splitEachLine("\t"){ toks->
    isoSet << toks[0]
}

reader.close()

println "${isoSet.size()} isoforms read."

def writer = new PrintWriter(outFile)
writer.println("snp\tisoform\tcorr\tpval\ttype")

def written = [] as TreeSet
int count = 0

new File(dir).eachDir{ fDir->
    def file = new File(fDir.absolutePath+'/'+CORREL_FILE)
    
    if( file.exists() ){
        count++
        reader = (file.name.endsWith('.gz')) ? 
            new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
            
        reader.readLine() //skip header
        
        reader.splitEachLine("\t"){ toks->
            def key = buildKey(toks[0], toks[1])
 
            if( (toks.size()>6) && (toks[1] in isoSet) && !(key in written) ){
                double corr = 0.0d
                double pval = 1.0d
                try{
                    corr = toks[2] as Double//corr
                    pval = toks[5] as Double//corrected pval
                } catch(e){
                    corr = 0.0d
                    pval = 1.0d
                }
                
                writer.println("${toks[0]}\t${toks[1]}\t${corr}\t${pval}\t${type}")
                written << key
            }
        }
        
        reader.close()
    }
}

writer.close()

println "${count} files processed."
println "${written.size()} eQTLs written."
println "Gzip out file ${outFile}."
"gzip -f ${outFile}".execute().waitFor()

println "End time: ${new Date()}"

return 0