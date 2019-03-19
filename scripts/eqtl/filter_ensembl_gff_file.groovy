#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

def ID_REG = /(ID=)(\w+):(\w+)/
def PARENT_REG = /(Parent=)(\w+):(\w+)/
def REG_LIST = [ID_REG, PARENT_REG]
def REG_F_ID = 3

def SEQ_NAMES = ((1..22).collect{it as String} + ['X','Y']) as TreeSet

def createReader = { file->
    (file.name.endsWith('.gz')) ? 
        new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
}

if( !this.args || this.args.length!=4 ){
    println 'Filter an ensembl GFF3 file by isoform or gene id.'
    println 'Usage: filter_ensembl_gff_file.groovy <gff_file> <filter_file> <column(s)> <out_file>'
    return 1
}

def inputGtf = this.args[0]
def filterFile = this.args[1]
def columns = (this.args[2].split(',') as List).collect{it as Integer}
def outFile = this.args[3]

println "Input : ${inputGtf}"
println "Output : ${outFile}"
println "Filter file: ${filterFile}"

//get list of terms
def reader = createReader(new File(filterFile))
def set = [] as TreeSet
reader.readLine()//skip header
reader.splitEachLine("\t"){ toks->
    columns.each{set << toks[it]} 
}

reader.close()

println "${set.size()} key terms read."


def file = new File(inputGtf)
reader = createReader(file)
        
def writer = new File(outFile).newWriter()
        
reader.eachLine{ line->
    if( line.startsWith('#') ){
        def toks = (line.split("\\s") as List).findAll{it}
        
        if( toks[0]=='##sequence-region' && (toks[1] in SEQ_NAMES) ){
            writer.println(line)
        }
        else if( toks[0]!='##sequence-region' ){
            writer.println(line)
        }
    }
    else{
        def toks = line.split("\t")
        //seqname source feature start end score strand frame attribute
        if( REG_LIST.any{ def mat=(toks[8]=~it);((mat) ? mat[0][REG_F_ID] in set : false) } ) {
            
            REG_LIST.each{ //clean Parent and ID attributes
                def mat = (line=~it)
                if(mat){
                    line = line.replace(mat[0][0], mat[0][1]+mat[0][REG_F_ID])
                }
            }
            
            if( toks[2]=='transcript' ){
                line = line.replace("\ttranscript\t","\tmRNA\t")
            }
            
            if( toks[0].startsWith('chr') ){
                writer.println(line)  
            }
            else{
                writer.println('chr'+line) 
            }
        }
    }
}

writer.close()
reader.close()

return 0