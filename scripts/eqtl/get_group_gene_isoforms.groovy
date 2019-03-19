#!/usr/bin/env groovy


//constants

String CORREL_UNFILT_FILE = 'correlation.txt'

//correlation.txt file fields
int F_CORR_SNP = 0
int F_CORR_ISO = 1
int F_CORR_CORR = 2
int F_CORR_PVAL = 5

//filtered isoforms file
int F_ISO_ID = 0
int F_ISO_GENE = 5

//best eqtls file
int F_BEST_ISO = 0
int F_BEST_SNP = 1

//end constants

def geneEqtls = [:] as TreeMap
def isoGeneMap = [:] as TreeMap

// closures

def geneSnpKey = { g, s-> "${g}:${s}" }

def splitKey = { key-> key.split(":") }

def doubleVal = { str-> try{return (str as Double)} catch(e){return 'NA'} }

def parseCorrFile = { file, chr->
    def reader = file.newReader()
    reader.readLine()//skip header

    reader.splitEachLine("\\s"){ toks->
        def gene = (toks.size()>5) ? isoGeneMap[toks[F_CORR_ISO]] : null
        if( gene ){
            def snp = toks[F_CORR_SNP]
            def key = geneSnpKey(gene, snp)

            if( geneEqtls.containsKey(key) ){
                geneEqtls[key] << [(F_CORR_ISO):toks[F_CORR_ISO],
                    (F_CORR_CORR):doubleVal(toks[F_CORR_CORR]),
                    (F_CORR_PVAL):doubleVal(toks[F_CORR_PVAL])]
            }
        }
    }

    reader.close()
}//

// end closures

if( !this.args || this.args.length!=4 ) {
    println "Groups isoforms by gene/snp eQTLs"
    println 'Usage: get_group_gene_isoforms.groovy <best eqtl filt. file> <isoform filt. file> <eqtl cis dir> <out file>'
    return 1
}

def bestEqtl = this.args[0]
def isoFile = this.args[1]
def cisDir = this.args[2]
def outFile = this.args[3]

println "Best eQTLs filtered file: ${bestEqtl}"
println "Isoforms filtered file: ${isoFile}"
println "eQTLs cis result dir: ${cisDir}"
println "Out file: ${outFile}"

println "Start time: ${new Date()}\n"

// 1 - read filtered isoforms file, fill iso->gene map
def reader = new File(isoFile).newReader()
reader.readLine()//skip header

reader.splitEachLine("\t"){ toks->
    isoGeneMap[toks[F_ISO_ID]] = toks[F_ISO_GENE]
}

reader.close()
println "${isoGeneMap.size()} isoforms read."

// 2 - read best eqtls file
reader = new File(bestEqtl).newReader()
reader.readLine()//skip header

reader.splitEachLine("\t"){ toks->
    def gene = isoGeneMap[toks[F_BEST_ISO]]
    def snp = toks[F_BEST_SNP]
    def key = geneSnpKey(gene, snp)
    
    if( !geneEqtls.containsKey(key) ){
        geneEqtls[key] = []
    }
}

reader.close()
println "${geneEqtls.size()} gene-snp eQTLs read."

// 3 - traverse eqtl cis results subfolders and print out result
new File(cisDir).eachDir{ dir->
    def file = new File(dir.absolutePath+'/'+CORREL_UNFILT_FILE)
    
    if( file.exists() ){
        def chr = dir.name.substring(0, dir.name.indexOf('_'))
        parseCorrFile(file, chr)
    }
}

def writer = new PrintWriter(outFile)
writer.println("gene\tsnp\tiso\tcorr\tpval")//header

geneEqtls.each{ key, list->
    def toks = splitKey(key)
    def written = [] as Set
    list.each{
        if( !(it[F_CORR_ISO] in written) ){
            writer.println("${toks[0]}\t${toks[1]}\t${it[F_CORR_ISO]}\t${it[F_CORR_CORR]}\t${it[F_CORR_PVAL]}")
            written << it[F_CORR_ISO]
        }
    }
}
writer.close()

println "\nEnd time: ${new Date()}\n"

return 0