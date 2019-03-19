#!/usr/bin/env groovy

/**
 * search isoform/snp that are eqtls cis and trans
 * 
 * NOT VALID, TODO: check whether the eqtl-trans satisfy trans condition
 */

String CORREL_FILE = 'correlation.txt.filter'

double pValThr = 1e-6
double pValThrTrans = 0.05

//correlation.txt file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 5

if( !this.args || this.args.length!=4 ){
    println 'Search eQTLs that are eqtls cis and trans'
    println 'Usage: search_eqtls_cis_trans.groovy <dir cis> <dir trans> (--snps|--isoforms) <out file>'
    return 1
}

String dirCis = this.args[0]
String dirTrans = this.args[1]
String mode = this.args[2]
String outFile = this.args[3]

println "input eqtl-cis dir: ${dirCis}"
println "input eqtl-trans dir: ${dirTrans}"
println "mode: ${mode}"
println "output file: ${outFile}"

def field = (mode=='--snps') ? F_SNP : F_ISO

// fill snps/isoforms list
def listEqtlCis = [] as TreeSet
def listEqtlCisTrans = [] as TreeSet

// closures

def fillEqtlCisList = { file->
    def reader = file.newReader()
    reader.readLine()

    reader.splitEachLine("\\s"){ toks->
        double pval = toks[F_PVAL] as Double
        
        if( pval<pValThr ) {
            listEqtlCis << toks[field]
        }
    }

    reader.close()
}

def fillEqtlCisTransList = { file->
    def reader = file.newReader()
    reader.readLine()

    reader.splitEachLine("\\s"){ toks->
        double pval = toks[F_PVAL] as Double
        
        if( (toks[field] in listEqtlCis) && pval<pValThrTrans ){
            listEqtlCisTrans << toks[field]
        }
    }

    reader.close()
}

// end closures

println "Start time: ${new Date()}\n"

// traverse eqtl-cis results dir
new File(dirCis).eachDir{ fDir->
    def file = new File(fDir.absolutePath+'/'+CORREL_FILE)
    if( file.exists() ) {
        fillEqtlCisList(file)
    }
}

println "${listEqtlCis.size()} eqtls-cis below pvalue ${pValThr}. ${new Date()}"

// traverse eqtl-trans results dir
new File(dirTrans).eachDir{ fDir->
    def file = new File(fDir.absolutePath+'/'+CORREL_FILE)
    if( file.exists() ){
        fillEqtlCisTransList(file)
    }
}

println "${listEqtlCisTrans.size()} eqtls trans/cis below pvalue ${pValThrTrans}. ${new Date()}"

//write output
def writer = new PrintWriter(outFile)
listEqtlCisTrans.each{ writer.println(it) }
writer.close()

println "End time: ${new Date()}"

return 0