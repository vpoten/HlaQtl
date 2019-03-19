#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

//constants
int F_SNP = 0
int F_ISO = 2
int F_LD = 9
def MISSING = ['null','NA']
int MIN_SIZE = 4


def createReader = { file->
    (file.name.endsWith('.gz')) ? 
        new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
}//

if( !this.args || this.args.length!=3 ){
    println 'Group missense eqtls by snps, discards eqtls with LD rSquare > threshold'
    println 'Usage: get_group_eqtl_missense.groovy <eqtl ld file> <ld_threshold> <outFile>'
    return 1
}

String inputFile = this.args[0]
Double ldThr = this.args[1] as Double
String outFile = this.args[2]

println "Input: ${inputFile}"
println "Output: ${outFile}"
println "Threshold: ${ldThr}"


def reader = createReader(new File(inputFile))
reader.readLine() //skip header

def groups = [:] as TreeMap
int rejected = 0
        
reader.splitEachLine("\t"){ toks->
    //fields: snp snp_locus iso iso_locus biotype corr pval snp_cis snp_cis_locus LD_rSq type
    boolean keep = true
    
    if( !(toks[F_LD] in MISSING) ){
        double rSq = toks[F_LD] as Double
        if( rSq>ldThr ){ keep = false; rejected++ }
    }
    
    if( keep ) {
        def snp = toks[F_SNP]
        def list = groups[snp]
        
        if( list==null ){
            list = []
            groups[snp] = list
        }
        
        list << [(F_ISO):toks[F_ISO]]
    }
}
        
reader.close()

int valid = groups.keySet().sum{ (groups[it].size()>MIN_SIZE) ? 1 : 0 }

println "${groups.size()} groups"
println "${valid} valid groups"
println "${rejected} rejected eqtls"

def writer = new PrintWriter(outFile)
//writer.println("isoform\tsnp\t\t\t") //write header

groups.each{ snp, list->
    if( list.size()>MIN_SIZE ){
        list.each{ writer.println("${it[F_ISO]}\t${snp}") }
    }
}

writer.close()

return 0