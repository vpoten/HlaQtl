#!/usr/bin/env groovy

def PLINK = 'plink'
def MISS_GENO = '--missing-genotype N'

//snp data map fields
def SNP_RISK_AL = 0
def SNP_OR = 1
def SNP_LOG_OR = 2

String ACT_1AL_RISK = '1risk'
String ACT_2AL_RISK = '2risk'
String ACT_ADDITIVE = 'additive'

def LIST_ACTION = [ ACT_1AL_RISK, ACT_2AL_RISK, ACT_ADDITIVE ]

if( !this.args || this.args.length!=4 ){
    println 'Transform a pair of ped/map files into an arff (weka) file'
    println 'usage: ped_to_arff <action> <ped_file> <map_file> <out_file>'
    println "Available actions: ${LIST_ACTION}"
    return 1
}

String action = this.args[0]

if( !(action in LIST_ACTION) ){
    println "Error: wrong action (${action})"
    return 1
}

String pedFile = this.args[1]
String mapFile = this.args[2]
String outFile = this.args[3]
    
    
println "input files: ${pedFile} / ${mapFile}"
println "output file: ${outFile}"
println "action: ${action}"

//// helper closures

def plinkTest = { test, ped, map, output->
    "${PLINK} --noweb --ped ${ped} --map ${map} ${test} --out ${output}"
}//


def getClearTokens = { string->
    (string.split("\\s") as List).findAll{!it.isEmpty()}
}//


def parseAssocResult = { file, map->
    def reader = new File(file).newReader()
    
    //get header and field indexes
    def header = reader.readLine().trim()
    def fields = getClearTokens(header)
    int snpIdx = fields.findIndexOf{it=='SNP'}
    int orIdx = fields.findIndexOf{it=='OR'}
    int a1Idx = fields.findIndexOf{it=='A1'}
    int a2Idx = fields.findIndexOf{it=='A2'}
         
    reader.eachLine{ line ->
        def toks = getClearTokens( line.trim() )
        
        def snpData = [(SNP_RISK_AL):toks[a1Idx]]

        try{
            double or = toks[orIdx] as Double
            
            if( or<1.0 ) {
                or = 1.0/or
                snpData[SNP_RISK_AL] = toks[a2Idx]
            }
            
            snpData[SNP_OR] = or
            snpData[SNP_LOG_OR] = Math.log(or)
            
            map[toks[snpIdx]] = snpData
            
        } catch(NumberFormatException e){}
    }
    
    reader.close()
}//


def calcRiskAllele = { snpData, al1, al2->
    def risk = snpData[SNP_RISK_AL]
    int sum = ((al1==risk) ? 1 : 0) + ((al2==risk) ? 1 : 0)
    
    if( action==ACT_1AL_RISK ){
        return (sum>0) ? 1 : 0
    }
    else if( action==ACT_2AL_RISK ){
        return (sum>1) ? 1 : 0
    }
    else if( action==ACT_ADDITIVE ){
        return sum
    }
        
    return 0
}//


def actionAttRange = {
    if( action==ACT_ADDITIVE ){
        return '{0,1,2}'
    }
    
    return '{0,1}'
}//

//// end of closures section

//read map file
def snpList = []

new File(mapFile).splitEachLine("\\s"){ toks->
    snpList << toks[1] // stores snp_id
}

//perform plink --assoc test
def comm = plinkTest("--assoc ${MISS_GENO}", pedFile, mapFile, outFile)

if( comm.execute().waitFor()!=0 ) {
    println "Error executing ${comm}"
    return 1
}

//parse test results, get risk allele
def snpMap = [:] as TreeMap
parseAssocResult(outFile+'.assoc', snpMap)

//generate arff file
def writer = new PrintWriter(outFile)

def fped = new File(pedFile)
def relation = fped.name

writer.println("@RELATION '${relation}'")
writer.println('')

snpList.each{ writer.println("@ATTRIBUTE '${it}' ${actionAttRange()}") }
writer.println("@ATTRIBUTE phenotype {1,2}")

writer.println('')
writer.println('@DATA')

fped.splitEachLine("\\s"){ toks->
    //PED-> Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype
    writer.println("%${toks[0]}:${toks[1]}")
    
    for(int i=6; i<toks.size(); i+=2) { 
        def snpId = snpList[(int)((i-6)/2)]
        def snpData = snpMap[snpId]
        writer.print("${calcRiskAllele(snpData,toks[i],toks[i+1])},") 
    }
    
    writer.println(toks[5])// print phenotype (class attribute)
}

writer.close()

return 0