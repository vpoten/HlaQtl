#!/usr/bin/env groovy

def MISSENSE_KEYS = ['41','42','44']
def CHR_CODES = ((1..22).collect{it as String} + ['X','Y']) as TreeSet
def FUN_FIELD = 11
def FREQ_FIELD = 3
def CHR_FIELD = 1
def SUFF_MISS = '.missense'
def SUFF_FILTER = '.freq_filter'
double freqThr = 0.01
double freqThrMax = 1.0-freqThr

if( !this.args || this.args.length!=4 ){
    println 'Filter a SNPContigLocusId.bcp from dbSNP to extract missense snps'
    println 'Usage: filter_missense_snps.groovy <bcp SNPContigLocusId> <bcp SNPAlleleFreq> <bcp SNPChrPosOnRef> <out file pref>'
    return 1
}

def file = this.args[0]
def freqFile = this.args[1]
def chrPosFile = this.args[2]
def outfile = this.args[3]

println "Input bcp file: ${file}"
println "Freqs. bcp file: ${freqFile}"
println "Output file: ${outfile}"
println "Freq. threshold: ${freqThr}"

println "Start time: ${new Date()}\n"

// 1- get missense SNPs
def reader = new BufferedReader( new FileReader(file) )
def writer = new PrintWriter(outfile+SUFF_MISS)

def missensSet = [] as TreeSet
        
reader.splitEachLine("\t"){ toks->
    if( toks[FUN_FIELD] in MISSENSE_KEYS ){
        Long id = toks[0] as Long
        if( !(id in missensSet) ) {
            writer.println( "rs${id}\t${toks[FUN_FIELD]}" )
            missensSet << id
        }
    }
}

reader.close()
writer.close()

println "#${missensSet.size()} missense snps"

// 2- filter by MAF freq
def filteredSet = [] as TreeSet
reader = new BufferedReader( new FileReader(freqFile) )

reader.splitEachLine("\t"){ toks->
    long id = toks[0] as Long
    double freq = toks[FREQ_FIELD] as Double
    
    if( (id in missensSet) && freq>freqThr && freq<freqThrMax )
        filteredSet << id
}

reader.close()

println "#${filteredSet.size()} missense snps after freq. filter"

// 3- Add chromosome
def chrSnpSet = [] as TreeSet
reader = new BufferedReader( new FileReader(chrPosFile) )

reader.splitEachLine("\t"){ toks->
    Long id = toks[0] as Long
    String chrId = toks[CHR_FIELD]
    
    if( (id in filteredSet) && (chrId in CHR_CODES) ){
        chrSnpSet << "chr${chrId}\trs${id}"
    }
}

reader.close()

println "#${chrSnpSet.size()} snps after chromosome filter"

// write filtered snps
writer = new PrintWriter(outfile+SUFF_FILTER)
chrSnpSet.each{ writer.println(it) }
writer.close()

println "End time: ${new Date()}\n"

return 0
  