#!/usr/bin/env groovy

def MISSING_CODES = ['','NA','NR','-','null']
    
def snpsMap = [:] as TreeMap //key=snp, value=set of features
def featuresMap = [:] as TreeMap //key=feature, value=set of snps
    
/**
 *
 */
def splitField = { str->
    if( !str || str in MISSING_CODES ){
        return null
    }
    return (str.split(",") as List).collect{it.trim()}
}

/**
 *
 */
def addToMap = {id, map, list->
    if( !map.containsKey(id) ){
        map[id] = [] as Set
    }
    map[id].addAll(list)
}

/**
 *
 */
def addFeatures = { chr, snpsF, featsF->
    def listSnps = splitField(snpsF)
    def listFeats = splitField(featsF)

    if( listSnps && listFeats ){
        listSnps = listSnps.findAll{it.startsWith('rs')}.collect{ it.contains(':') ? it.split(':') as List : it }.flatten()
        listSnps.each{snp-> addToMap(snp,snpsMap,listFeats) }
        listFeats.each{feat-> addToMap(feat,featuresMap,listSnps) }
    }
}

if( !this.args || this.args.length!=2 ){
    println 'Reads a gwas catalog file and converts it to snp - features format'
    println 'Usage: gwas_catalog_filter.groovy <gwas-catalog file> <out_file>'
    return 1
}

def fileIn = this.args[0]
def fileOut = this.args[1]

def reader = new File(fileIn).newReader()
def header = reader.readLine() //skip header
def fields = header.split("\t") as List

def featIdx = fields.findIndexOf{it=='Disease/Trait'}
def chrIdx = fields.findIndexOf{it=='Chr_id'}
def snpIdx = fields.findIndexOf{it=='SNPs'}


reader.splitEachLine("\t"){ toks->
    def chr = toks[chrIdx]
    def snps = toks[snpIdx]
    def feats = toks[featIdx]
    addFeatures(chr, snps, feats)
}

reader.close()

def writer = new PrintWriter(fileOut)
writer.println("snp\tfeatures")//write output file header

snpsMap.each{ snpId, feats->
    writer.println("${snpId}\t${feats.sum{it+','}}\t")
}

writer.close()

return 0