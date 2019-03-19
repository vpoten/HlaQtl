#!/usr/bin/env groovy

//correlation filtered file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 3
int F_TYPE = 4

def corrRates = [2*0.012, 2*0.0085, 2*0.007, 2*0.006, 2*0.006, 2*0.006]
def snpsMaf = [:] as TreeMap //snps frequencies map

// closures
def buildEqtlData = { toks->
    return [ (F_SNP):toks[F_SNP], 
        (F_CORR):Math.abs(toks[F_CORR] as Double),
        (F_PVAL):(toks[F_PVAL] as Double),
        (F_TYPE):toks[F_TYPE] ]
}

def calcDiff = { snpId, corr->
    double maf = snpsMaf[snpId]
    int idx = maf/0.1
    return corrRates[idx]
}

// end closures


if( !this.args || this.args.length!=3 ) {
    println "Get group of best eqtls from a eqtl all filtered file, threshold depends on MAF"
    println 'Usage: get_group_best_eqtl_maf.groovy <eqtl all filtered file> <snps info log> <out file>'
    return 1
}

def input = this.args[0]
def snpLog = this.args[1]
def outfile = this.args[2]

println "Input file: ${input}"
println "SNPs log file: ${snpLog}"
println "Out file: ${outfile}"
println "Correlation similarity rates: ${corrRates}"

println "Start time: ${new Date()}\n"

// parse snps log file
reader = new File(snpLog).newReader()
reader.readLine()//skip header

reader.splitEachLine("\t"){ toks->
    snpsMaf[toks[0]] = toks[6] as Double
}

reader.close()

println "${snpsMaf.size()} SNPs freq. read."


//get best eQTLs
reader = new File(input).newReader()
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
        double diff = (best) ? calcDiff(best[F_SNP], best[F_CORR]) : 0.0
        
        if( best && ((Math.abs(best[F_CORR]-corr) < diff)) ){
            groupBest[toks[F_ISO]] <<  buildEqtlData(toks)
        }
    }
}

reader.close()

// write groups of best eqtls
def writer = new PrintWriter(outfile)
writer.println("isoform\tsnp\tcorrelation\tpval\ttype\tmaf")//write header

groupBest.each{ iso, list->
    def best = bestEqtls[iso]
    writer.println("${iso}\t${best[F_SNP]}\t${best[F_CORR]}\t${best[F_PVAL]}\t${best[F_TYPE]}\t${snpsMaf[best[F_SNP]]}")
    
    list.each{
        writer.println("${iso}\t${it[F_SNP]}\t${it[F_CORR]}\t${it[F_PVAL]}\t${it[F_TYPE]}\t${snpsMaf[it[F_SNP]]}")
    }
}

writer.close()

println "\nEnd time: ${new Date()}\n"

return 0
