#!/usr/bin/env groovy

// group best eqtl file fields
int GBE_ISO = 0
int GBE_SNP = 1
int GBE_CORR = 2
int GBE_PVAL = 3
int GBE_TYPE = 4

// group gene-snp file fields
int GG_GENE = 0
int GG_SNP = 1
int GG_ISO = 2
int GG_CORR = 3
int GG_PVAL = 4

def buildEqtlDataGG = { toks->
    return [ (GG_CORR):Math.abs(toks[GG_CORR] as Double),
        (GG_PVAL):(toks[GG_PVAL] as Double), (GG_ISO):toks[GG_ISO] ]
}

if( !this.args || this.args.length!=3 ) {
    println "Join gene-snp eqtl groups where snps are best-equivalent for any gene transcript of group."
    println 'Usage: join_gene_snp_groups.groovy <group gene snp file> <group best eqtls file> <out file>'
    return 1
}


def groupGeneFile = this.args[0]
def groupBestFile = this.args[1]
def outFile = this.args[2]

println "Input group gene-snp file: ${groupGeneFile}"
println "Input group eqtl best file: ${groupBestFile}"
println "Out file: ${outFile}"

println "Start time: ${new Date()}\n"

// 1 - get best eQTLs groups
def reader = new File(groupBestFile).newReader()
reader.readLine()//skip header

def bestEqtls = [:] as TreeMap// map with key=isoId, value=set of snps

reader.splitEachLine("\t"){ toks->
    def iso = toks[GBE_ISO]
    def list = bestEqtls[iso]
    if( list==null ){
        list = [] as TreeSet
        bestEqtls[iso] = list
    }
    
    list << toks[GBE_SNP]
}

reader.close()

println "${bestEqtls.size()} isoforms read"

// 2 - read gene eqtls groups file
def geneEqtls = [:] as TreeMap

reader = new File(groupGeneFile).newReader()
def header = reader.readLine()//skip header

reader.splitEachLine("\t"){ toks->
    def gene = toks[GG_GENE]
    def map = geneEqtls[gene]
    if( map==null ){
        map = [:]
        geneEqtls[gene] = map
    }
    
    def list = map[toks[GG_SNP]]
    if( list==null ){
        list = []
        map[toks[GG_SNP]] = list
    }
    
    list << buildEqtlDataGG(toks)
}

reader.close()

println "${geneEqtls.size()} genes read"


// 3 - join
geneEqtls.each{ gene, map->
    if( map.size()>1 ){
        def snps = map.keySet() as List
        def isoforms = map[snps[0]].collect{ it[GG_ISO] }
        
        // get subsets of equivalent snps
        def joins = isoforms.collect{iso-> snps.findAll{it in bestEqtls[iso]} as Set}
        joins = joins.findAll{it?.size()>1}
        
        // search for transitive relations and join subsets
        boolean joined = true
        
        while( joined ){
            joined = false
            
            for(i=1; i<joins.size(); i++){
                if( joins[i-1].any{it in joins[i]} ){
                    println "joined ${joins[i-1]} and ${joins[i]}"
                    joins[i-1]+=joins[i]
                    joins.remove(i)
                    joined = true
                    break
                }
            }
        }
        
        // add joined subsets of snps and remove single snp eqtls
        joins.each{ snpSet->
            def newSnp = snpSet.sum{it+','}
            newSnp = newSnp.substring(0,newSnp.length()-1)
            
            map[newSnp] = isoforms.collect{iso-> snpSet.collect{snp-> map[snp].find{it[GG_ISO]==iso}}.max{it[GG_CORR]}}
            snpSet.each{ map.remove(it) }
        }
    }
}

// 4 - write output
def writer = new PrintWriter(outFile)
writer.println(header)//write header

geneEqtls.each{ gene, map->
    map.each{ snp, list->
        list.each{ 
            writer.println("${gene}\t${snp}\t${it[GG_ISO]}\t${it[GG_CORR]}\t${it[GG_PVAL]}")
        }
    }
}

writer.close()


println "\nEnd time: ${new Date()}\n"

return 0