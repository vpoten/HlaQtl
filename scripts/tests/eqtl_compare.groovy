#!/usr/bin/env groovy

//correlation.txt file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 5

if( !this.args || this.args.length!=2 ){
    println 'Eqtl result file comparison'
    println 'Usage: eqtl_compare.groovy <file1> <file2>'
    return 1
}

//// closures

def buildEqtlData = { toks->
    return [ (F_SNP):toks[F_SNP], 
        (F_CORR):Math.abs(toks[F_CORR] as Double),
        (F_PVAL):(toks[F_PVAL] as Double),
        (F_ISO):toks[F_ISO] ]
}
// end closure

def buildKey = { snp, iso->
    "${snp}:${iso}"
}// end closure

////

def files = [this.args[0], this.args[1]]
def mapEqtls = [([:] as TreeMap), ([:] as TreeMap)]

def snps = [] as TreeSet
def isoforms = [] as TreeSet


(0..1).each{ i->
    def reader = new File(files[i]).newReader()
    reader.readLine()//skip header
    def map = mapEqtls[i]

    reader.splitEachLine("\t"){ toks->
            snps << toks[F_SNP]
            isoforms << toks[F_ISO]
            map[buildKey(toks[F_SNP],toks[F_ISO])] = buildEqtlData(toks)
    }

    reader.close()
}

def writer = System.out

snps.each{ snp->
    isoforms.each{ iso->
        def key = buildKey(snp, iso)
        def eqtl1 = mapEqtls[0][key]
        def eqtl2 = mapEqtls[1][key]
        
        if( eqtl1 || eqtl2) {
            writer.println("${snp}\t${iso}\t${eqtl1?.get(F_CORR)}\t${eqtl1?.get(F_PVAL)}\t${eqtl2?.get(F_CORR)}\t${eqtl2?.get(F_PVAL)}")
        }
    }
}

return 0
