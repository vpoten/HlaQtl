#!/usr/bin/env groovy

def CHR_CODES = ((1..22).collect{it as String} + ['X','Y']) as TreeSet
def CHR_FIELD = 1
def POS_FIELD = 2

if( !this.args || this.args.length!=3 ){
    println 'Uses a SNPChrPosOnRef.bcp from dbSNP to add rs ids to genome positions'
    println 'Usage: add_snpid.groovy <positions file> <bcp SNPChrPosOnRef> <out file>'
    return 1
}

def file = this.args[0]
def chrPosFile = this.args[1]
def outfile = this.args[2]

println "Input chr position: ${file}"
println "Snp position. bcp file: ${chrPosFile}"
println "Output file: ${outfile}"

println "Start time: ${new Date()}\n"

//parse chr position file
def positions = [:] as TreeMap

def reader = new File(file).newReader()
reader.readLine() //skip header

reader.splitEachLine("\t"){ toks->
    def chr = toks[0].startsWith('chr') ? toks[0].substring(3) : toks[0]
    def position = toks[1] as Integer
    
    def set = positions[chr]
    
    if(set==null){
        set = [] as TreeSet
        positions[chr] = set
    }
    
    set << position
}

reader.close()


// add rs id comparing using chromosome and position
def writer = new PrintWriter(outfile)
writer.println("snp\tchr\tposition")//write output file header

def chrSnpSet = [] as TreeSet
reader = new File(chrPosFile).newReader()

reader.splitEachLine("\t"){ toks->
    Long id = toks[0] as Long
    String chrId = toks[CHR_FIELD]
    
    if( toks[POS_FIELD] && chrId in CHR_CODES ){
        Integer pos = toks[POS_FIELD] as Integer

        if( positions.containsKey(chrId) && (pos in positions[chrId]) ){
            writer.println("rs${id}\t${chrId}\t${pos}")
        }
    }
}

reader.close()
writer.close()

println "End time: ${new Date()}\n"

return 0