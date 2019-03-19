#!/usr/bin/env groovy

//file fields

// group best eqtl
int GBE_ISO = 0
int GBE_SNP = 1
int GBE_CORR = 2
int GBE_PVAL = 3
int GBE_TYPE = 4

//eqtl all
int EA_SNP = 0
int EA_ISO = 1
int EA_CORR = 2
int EA_PVAL = 3
int EA_TYPE = 4

// eqtl info
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 3
int F_TYPE = 4

//group labels
def LBL_G1 = 'best1'
def LBL_G2 = 'best2'
def LBL_G3 = 'best3'
def LBL_NB1 = 'notbest1'
def LBL_NBW = 'notbestw'

//alternative
def LBL_BEST = 'best'
def LBL_NOTBEST = 'notbest'

boolean binaryGroups = true //if true groups are classified int best and notbest
boolean removeUniq = false //set true for removing of unique best eqtls

// eqtl comparator (desc)
def eqtlComp = [compare:{a,b-> b[(F_CORR)]<=>a[(F_CORR)]}] as Comparator

// closures
def buildEqtlDataGBE = { toks->
    return [ (F_SNP):toks[GBE_SNP], (F_CORR):Math.abs(toks[GBE_CORR] as Double),
        (F_PVAL):(toks[GBE_PVAL] as Double), (F_TYPE):toks[GBE_TYPE] ]
}

def buildEqtlDataEA = { toks->
    return [ (F_SNP):toks[EA_SNP], (F_CORR):Math.abs(toks[EA_CORR] as Double),
        (F_PVAL):(toks[EA_PVAL] as Double), (F_TYPE):toks[EA_TYPE] ]
}

// end closures

def snpGroups

if( binaryGroups ){
    snpGroups = [(LBL_BEST):[] as TreeSet, (LBL_NOTBEST):[] as TreeSet]
}
else{
    snpGroups = [(LBL_G1):[] as TreeSet, (LBL_G2):[] as TreeSet, (LBL_G3):[] as TreeSet, 
        (LBL_NB1):[] as TreeSet, (LBL_NBW):[] as TreeSet]
}


if( !this.args || this.args.length!=4 ) {
    println "Generates an extended snps info log file with new attributes of importance in eqtl groups."
    println 'Usage: get_snps_groups.groovy <group best eqtls file> <eqtl all filtered file> <snps info log> <out file>'
    return 1
}

def groupBestFile = this.args[0]
def eqtlAllFile = this.args[1]
def snpLog = this.args[2]
def outFile = this.args[3]

println "Input group eqtl best file: ${groupBestFile}"
println "eQTL all filtered file: ${eqtlAllFile}"
println "SNPs log file: ${snpLog}"
println "Out file: ${outFile}"

println "Start time: ${new Date()}\n"

// 1 - get best eQTLs groups
def reader = new File(groupBestFile).newReader()
reader.readLine()//skip header

def bestEqtls = [:] as TreeMap// map with key=isoId, value=list of eqtlInfo

reader.splitEachLine("\t"){ toks->
    def iso = toks[GBE_ISO]
    def list = bestEqtls[iso]
    if( list==null ){
        list = []
        bestEqtls[iso] = list
    }
    
    list << buildEqtlDataGBE(toks)
}

reader.close()
println "${bestEqtls.size()} isoforms read in best eQTLs"

if( removeUniq ){
    def uniq = bestEqtls.findAll{ it.value.size()==1 }
    uniq.each{ bestEqtls.remove(it.key) }
    println "${uniq.size()} unique best eqtls removed"
}

//sort groups of eqtls by correlation (descend) and add to groups
bestEqtls.each{ iso, list-> 
    Collections.sort(list,eqtlComp)
    if( binaryGroups ){
        list.each{ snpGroups[(LBL_BEST)] << it[F_SNP] }
    }
    else{
        snpGroups[(LBL_G1)] << list.first()[F_SNP]
        if( list.size()>1 ){ snpGroups[(LBL_G2)] << list[1][F_SNP] }
        if( list.size()>2 ){ snpGroups[(LBL_G3)] << list[2][F_SNP] }
    }
}



// 2 - get eQTLs not present in best groups
reader = new File(eqtlAllFile).newReader()
reader.readLine()//skip header

def notBestEqtls = [:] as TreeMap// map with key=isoId, value=list of eqtlInfo

reader.splitEachLine("\t"){ toks->
    def iso = toks[EA_ISO]
    def type = toks[EA_TYPE]
    def listBest = bestEqtls[iso]
    double corr = Math.abs(toks[EA_CORR] as Double)
    
    if( listBest!=null && type!='trans' && listBest.last()[F_CORR]>corr ){
        def list = notBestEqtls[iso]
        if( list==null ){
            list = []
            notBestEqtls[iso] = list
        }
    
        list << buildEqtlDataEA(toks)
        
        if( binaryGroups ){
            //print out notbest eqtl
            println "${toks[0]}\t${toks[1]}\t${toks[2]}\t${toks[3]}\t${toks[4]}\tnotbest"
        }
    }
}

reader.close()
println "${notBestEqtls.size()} isoforms read in not best eQTLs"

//sort groups of eqtls by correlation (descend) and add to groups
notBestEqtls.each{ iso, list-> 
    Collections.sort(list,eqtlComp)
    if( binaryGroups ){
        list.each{ snpGroups[(LBL_NOTBEST)] << it[F_SNP] }
    }
    else{
        snpGroups[(LBL_NB1)] << list.first()[F_SNP]
        if( list.size()>1 ){ snpGroups[(LBL_NBW)] << list.last()[F_SNP] }
    }
}



// 3 - write output
reader = new File(snpLog).newReader()
def header = reader.readLine()//skip header

def writer = new PrintWriter(outFile)
writer.println("${header}\tattributes")//write header

reader.splitEachLine("\t"){ toks->
    def snp = toks[0]
    def atts = [] as TreeSet
    
    if( toks[1]=='23' ){ toks[1]='X' }
   
    snpGroups.each{ lbl, set-> 
        if( snp in set ){ atts << lbl }
    }
    
    if( atts ){
        (0..6).each{ writer.print("${toks[it]}\t") }
        writer.println( atts.sum{it+','} ) 
    }
}

writer.close()
reader.close()

println "\nEnd time: ${new Date()}\n"

return 0