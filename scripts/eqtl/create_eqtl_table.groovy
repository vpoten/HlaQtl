#!/usr/bin/env groovy

/**
 * create eqtl table file for graph generation
 */

def buildKey = { snp, iso->
    "${snp}_${iso}"
}//

def outFile = '/home/victor/Escritorio/eqtls_snps/all_eqtl_cis_trans.txt'
def allCis = '/home/victor/Escritorio/eqtls_snps/eqtl_cis.all.txt'
def bestCis = '/home/victor/Escritorio/eqtls_snps/best_eqtls_1e05.raw.out'
def transEqtl = '/home/victor/Escritorio/eqtls_snps/correlation.txt.trans.filter'

def writer = new PrintWriter(outFile)

writer.println("snp\tisoform\tcorrelation\tpvalue\ttype")

// 1 - read best eqtl cis
def reader = (new File(bestCis)).newReader()
def bestEqtls = [] as TreeSet

reader.splitEachLine("\\s"){ toks->
    bestEqtls << buildKey(toks[1],toks[0])
}

reader.close()

// 2 - all eqtl cis
reader = (new File(allCis)).newReader()

reader.splitEachLine("\\s"){ toks->
    def key = buildKey(toks[0],toks[1])
    def type = (key in bestEqtls) ? 'best' : 'cis'
    //snp\tisoform\tcorrelation\tpvalue\ttype
    writer.println("${toks[0]}\t${toks[1]}\t${toks[2]}\t${toks[5]}\t${type}")
}

reader.close()

// 3 - eqtl trans
reader = (new File(transEqtl)).newReader()
reader.readLine()//skip header

reader.splitEachLine("\\s"){ toks->
    //snp\tisoform\tcorrelation\tpvalue\ttype
    writer.println("${toks[0]}\t${toks[1]}\t${toks[2]}\t${toks[5]}\ttrans")
}

reader.close()

writer.close()

return 0