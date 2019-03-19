#!/usr/bin/env groovy

//correlation filtered file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 3
int F_TYPE = 4

double corrRate = 1e-6

// closure
def buildEqtlData = { toks->
    return [ (F_SNP):toks[F_SNP], 
        (F_CORR):Math.abs(toks[F_CORR] as Double),
        (F_PVAL):(toks[F_PVAL] as Double),
        (F_TYPE):toks[F_TYPE] ]
}
// end closure


if( !this.args || this.args.length!=2 ) {
    println "Get group of best eqtls from a eqtl all filtered file"
    println 'Usage: get_group_best_eqtl.groovy <eqtl all filtered file> <out file>'
    return 1
}

def input = this.args[0]
def outfile = this.args[1]

println "Input file: ${input}"
println "Out file: ${outfile}"
println "Correlation similarity rate: ${corrRate}"

println "Start time: ${new Date()}\n"

def snpsSet = [] as TreeSet

//get best eQTLs
def reader = new File(input).newReader()
reader.readLine()//skip header

def bestEqtls = [:] as TreeMap

reader.splitEachLine("\t"){ toks->
    if( toks[F_TYPE]=='best' ) {
        bestEqtls[toks[F_ISO]] = buildEqtlData(toks)
    }
}

reader.close()

println "${bestEqtls.size()} best eQTLs read."

//search group of best eQTLs among cis eQTLs
def groupBest = [:] as TreeMap
bestEqtls.each{ iso, v-> groupBest[iso] = [] }

reader = new File(input).newReader()
reader.readLine()//skip header

reader.splitEachLine("\t"){ toks->
    if( toks[F_TYPE]=='cis' ) {
        def best = bestEqtls[toks[F_ISO]]
        double corr = Math.abs(toks[F_CORR] as Double)
        
        if( best && (Math.abs(best[F_CORR]-corr) < best[F_CORR]*corrRate) ){
            groupBest[toks[F_ISO]] <<  buildEqtlData(toks)
        }
    }
}

reader.close()

// write groups of best eqtls
def writer = new PrintWriter(outfile)
writer.println("isoform\tsnp\tcorrelation\tpval\ttype")//write header

groupBest.each{ iso, list->
    def best = bestEqtls[iso]
    writer.println("${iso}\t${best[F_SNP]}\t${best[F_CORR]}\t${best[F_PVAL]}\t${best[F_TYPE]}")
    
    list.each{
        writer.println("${iso}\t${it[F_SNP]}\t${it[F_CORR]}\t${it[F_PVAL]}\t${it[F_TYPE]}")
    }
}

writer.close()

println "\nEnd time: ${new Date()}\n"

return 0