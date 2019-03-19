#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

def GENEID_REG = /gene_id \"([\w\.]+)\"/
def TRANSID_REG = /transcript_id \"([\w\.]+)\"/
def REG_LIST = [TRANSID_REG, GENEID_REG]

if( !this.args || this.args.length!=4 ){
    println 'Filter a GTF file by isoform or gene id.'
    println 'Usage: filter_gtf_file.groovy <gtf_file> <filter_file> <column> <out_file>'
    return 1
}

def inputGtf = this.args[0]
def filterFile = this.args[1]
int column = this.args[2] as Integer
def outFile = this.args[3]

println "Input : ${inputGtf}"
println "Output : ${outFile}"
println "Filter file: ${filterFile}"

//get list of terms
reader = new File(filterFile).newReader()
def set = [] as TreeSet
reader.readLine()//skip header
reader.splitEachLine("\t"){ toks->
    set << toks[column]
}

reader.close()

println "${set.size()} key terms read."


def file = new File(inputGtf)
def reader = (file.name.endsWith('.gz')) ? 
        new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
        
def writer = new File(outFile).newWriter()
        
reader.eachLine{ line->
    def toks = line.split("\t")
    //seqname source feature start end score strand frame attribute
    if( REG_LIST.any{ def mat=(toks[8]=~it);((mat) ? mat[0][1] in set : false) } ) {
          writer.println(line)  
    }
}

writer.close()
reader.close()

return 0