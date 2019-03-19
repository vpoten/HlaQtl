#!/usr/bin/env groovy

import org.ngsutils.Utils

//available actions
String ACT_LIST = 'list_regions'
String ACT_STATS = 'stats'
String ACT_REP_FAIL = 'repeat_failed'
String ACT_CALC = 'calc'
String ACT_LIST_FAIL = 'list_failed'
String ACT_GEN_FEAT = 'gen_features'
String ACT_SPLIT_FAIL = 'split_failed'
def LIST_ACTION = [ ACT_LIST, ACT_STATS, ACT_REP_FAIL, ACT_CALC, ACT_LIST_FAIL, 
    ACT_GEN_FEAT, ACT_SPLIT_FAIL ]

//constants
int THREADS = 10
int REG_STEP = 1000000
int REG_OVERLAP = 100000
String S3_VCF_FILE = 'ALL.{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz'
String CORREL_FILE = 'correlation.txt.filter'
String SNPS_FILE = '1000g_snps.map'

def GFF3_OUT = 'eqtls_tracks.gff3'
def BEDG_OUT = 'eqtls_tracks.bed'

//constants for region data map
def RD_REGION = 'region'
def RD_COMM = 'command'
def RD_DIR = 'dir'
def RD_OUTFILE = 'output'

//constants for snp eqtl info data map
def SNP_POS = 1
def SNP_ID = 2
def SNP_PVAL = 3
def SNP_CORR = 4
def SNP_ISO = 5
def SNP_EQTLS = 6

//correlation.txt file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 5

//hg19 chromosomes sizes
def chromSizes = [
'chr1':249250621,
'chr2':243199373,
'chr3':198022430,
'chr4':191154276,
'chr5':180915260,
'chr6':171115067,
'chr7':159138663,
'chr8':146364022,
'chr9':141213431,
'chr10':135534747,
'chr11':135006516,
'chr12':133851895,
'chr13':115169878,
'chr14':107349540,
'chr15':102531392,
'chr16':90354753,
'chr17':81195210,
'chr18':78077248,
'chr19':59128983,
'chr20':63025520,
'chr21':48129895,
'chr22':51304566,
'chrX':155270560,
///'chrY':59373566
]

if( !this.args || this.args.length!=1 ){
    println 'Genome eqtl calculation and statistics'
    println 'Usage: genome_eqtl.groovy <action>'
    println "Available actions: ${LIST_ACTION}"
    return 1
}

String action = this.args[0]

if( !(action in LIST_ACTION) ){
    println "Error: wrong action (${action})"
    return 1
}

// directories and options
String jarFile = '/home/users/ipb/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar'
String javaOpts = '-Xmx6g'
String awsProps = '/home/users/ipb/HlaQtl/ngsengine_fuen.properties'
String outDir = '/home/users/ipb/eqtl_out/'
String exprDataDir = '/home/users/ipb/cufflinks_ensOnly/'
String groups = 'CEU,GBR,FIN,TSI'
String cuffPre = 'cufflinks_ensOnly_'
String vcfDir = '/home/users/ipb/vcfFiles/'

println "Num. threads: ${THREADS}"
println "Step: ${REG_STEP}"
println "Overlap: ${REG_OVERLAP}"
println "Output dir: ${outDir}"
println "Expr. data dir: ${exprDataDir}"
println "Groups: ${groups}"

println "Start time: ${new Date()}\n"

int incr = REG_STEP - REG_OVERLAP
def regionsData = []

// createRegData closure 
def createRegData = { chr, start, end->
    def dir = "${chr}_${start}_${end}"
    def locus = "${chr}:${start}-${end}"
    def vcfFile = S3_VCF_FILE.replace('{chr}',chr)
     
    def command = "java ${javaOpts} -jar ${jarFile} eqtl --props=${awsProps} " +
        "--output=${outDir+dir} --locus=${locus} --group=${groups} " +
        "--cuffPre=${cuffPre} --exprData=${exprDataDir} --vcf=${vcfDir+vcfFile}"
        
    def outFile = (outDir+dir+'.out')
    
    return [ (RD_REGION):locus, (RD_COMM):command, 
        (RD_DIR):outDir+dir, (RD_OUTFILE):outFile ]
}// end of createRegData

// create regions to calculate
chromSizes.each{ chr, size->
    for(int start=0; start<size; start+=incr) {
        int end = ((start+REG_STEP)<size) ? start+REG_STEP : size
        
        def data = createRegData(chr,start,end)
        
        //mkdir output dir
        //"mkdir ${data[RD_DIR]}".execute().waitFor()
        regionsData << data
    }
}

//closure for failed tasks listing
def failedTasksArgs = { listRegs->
    def fails = []
    
    listRegs.each {
        if( !(new File(it[RD_OUTFILE]).exists()) ){
            fails << it
        }
        else{
            def proc = "tail -n 2 ${it[RD_OUTFILE]}".execute()
            proc.waitFor()
            if( !proc.text.contains('End time:') )
                fails << it //not finished ok
        }
    }
    
    return fails
}

def failedTasks = {
    return failedTasksArgs(regionsData)
}// end of failedTasks



//do actions
if( action==ACT_CALC ) {
    println "${regionsData.size()} eqtl jobs to execute:"
    Utils.runCommands(regionsData.collect{it[RD_COMM]}, THREADS, regionsData.collect{it[RD_OUTFILE]}, false)
}
else if( action==ACT_LIST ) {
    println "===== List of regions ====="
    regionsData.each{ println it[RD_REGION] }
    println "===== List of commands to execute ====="
    regionsData.each{ println it[RD_COMM] }
}
else if( action==ACT_STATS ) {
    def fails = failedTasks()
    println "Statistics:"
    println "Total ${regionsData.size()}"
    println "Failed ${fails.size()}"
    int corrReg = 0
    
    regionsData.each{
        def file = it[RD_DIR]+'/'+CORREL_FILE
        if( (new File(file).exists()) && Utils.countLines(file)>1 )
            corrReg++
    }
    
    println "Regions where corr. found ${corrReg}"
}
else if( action==ACT_REP_FAIL ) {
    def fails = failedTasks()
    println "${fails.size()} failed eqtl jobs to re-execute:"
    Collections.shuffle(fails, new Random())
    int threads = 6
    Utils.runCommands(fails.collect{it[RD_COMM]}, threads, fails.collect{it[RD_OUTFILE]}, false)
}
else if( action==ACT_LIST_FAIL ) {
    def fails = failedTasks()
    println "===== List of failed regions ====="
    fails.each{ println it[RD_REGION] }
    println "===== List of failed commands ====="
    fails.each{ println it[RD_COMM] }
}
else if( action==ACT_GEN_FEAT ) {
    //generate features tracks for eqtl results: bedGraph and gff3
    double pValThr = 0.001
    double corrThr = 0.1
    def currChr = null
    def chrSnpsMap = null
    
    def fails = failedTasks()
    def completeRegions = regionsData.findAll{ !(it in fails) }
    
    println "Generating track files for ${completeRegions.size()} regions"
    println "${fails.size()} failed eqtls"
    println "gff3 out: ${outDir+GFF3_OUT}"
    println "bedGraph out: ${outDir+BEDG_OUT}"
    
    //create writers for gff3 and bedGraph
    def writerGff = new File(outDir+GFF3_OUT).newPrintWriter()
    def writerBed = new File(outDir+BEDG_OUT).newPrintWriter()
    
    //write headers
    writerBed.println('track type=bedGraph name="eqtl_scores" description="eQTLs correlation values"')
    writerGff.println('##gff-version 3')
    
    def dumpChrData = { chr->
        //dump current chr data
        def list = chrSnpsMap.values() as List
        list = list.sort{ it[SNP_POS] }

        list.each{
            if( it[SNP_CORR]>corrThr ) {
                def attrs = "ID=${it[SNP_ID]};Name=${it[SNP_ID]};Note=eQTLs ${it[SNP_ISO]}%3D${it[SNP_CORR]} ${it[SNP_EQTLS].toString()}"
                writerGff.println("${chr}\teqtl\tSNP\t${it[SNP_POS]}\t${it[SNP_POS]}\t${it[SNP_CORR]}\t+\t.\t${attrs}")
                writerBed.println("${chr}\t${it[SNP_POS]}\t${it[SNP_POS]}\t${it[SNP_CORR]}")
            }
        }
    }//end closure
    
    completeRegions.eachWithIndex{ reg, i->
        def chr = (reg[RD_REGION]=~Utils.locusRegex)[0][1]
        println "Processing ${reg[RD_REGION]} ${i+1}/${completeRegions.size()}"
        
        if( currChr!=chr ){
            if( currChr!=null ){
                dumpChrData(currChr)
            }
            
            chrSnpsMap = [:] as TreeMap
            currChr = chr
        }
        
        // parse region .map file
        new File(reg[RD_DIR]+'/'+SNPS_FILE).eachLine{ line->
            def toks = line.split("\\s")
            if( !chrSnpsMap[toks[1]] )
                chrSnpsMap[toks[1]] = [ (SNP_ID):toks[1], (SNP_POS):toks[3] as Integer,
                    (SNP_CORR):0.0, (SNP_EQTLS):new StringBuilder() ]
        }
        
        //parse correlation file
        def reader = new File(reg[RD_DIR]+'/'+CORREL_FILE).newReader()
        reader.readLine()//skip header
        
        reader.splitEachLine("\\s"){ toks->
            double pval = toks[F_PVAL] as Double

            if( pval<pValThr ){
                //store snp-iso correlation value
                def snpData = chrSnpsMap[toks[F_SNP]]
                double corrVal = Math.abs( toks[F_CORR] as Double )
                
                if( snpData && snpData[SNP_CORR]<corrVal ){
                    snpData[SNP_CORR] = corrVal
                    snpData[SNP_PVAL] = pval
                    snpData[SNP_ISO] = toks[F_ISO]
                }
                if( snpData ){
                    snpData[SNP_EQTLS] << "%2C${toks[F_ISO]}%3D${toks[F_CORR]}"
                }
            }
        }
        
        reader.close()
    }
    
    dumpChrData(currChr)// dump data of the last chromosome
    
    [writerGff,writerBed].each{ it.close() }
    
}//end action ACT_GEN_FEAT
else if( action==ACT_SPLIT_FAIL ) {
    def fails = failedTasks()
    
    def failsSplit = [1:fails]
    def totals = [1:regionsData]
    
    [2,4].each{ parts->
        int step = REG_STEP/parts + REG_OVERLAP/parts
        totals[parts] = []
        
        failsSplit[(parts/2) as Integer].each{ reg->
            def mat = (reg[RD_REGION]=~Utils.locusRegex)
            def chr = mat[0][1]
            def start = mat[0][2] as Integer
            def end = mat[0][3] as Integer

            totals[parts] << createRegData(chr, start, start+step)
            totals[parts] << createRegData(chr, end-step, end)
        }
        
        failsSplit[parts] = failedTasksArgs(totals[parts])
    }
    
    //print summary
    [1,2,4].each{
        int step = REG_STEP/it + REG_OVERLAP/it
        println "Step: ${step}; ${totals[it].size()} tasks, ${failsSplit[it].size()} failed"
    }
    
    def failsToDo = failsSplit[4]
    
//    // create output dirs
//    failsToDo.each{ data->
//        if( !(new File(data[RD_DIR]).exists()) )
//            "mkdir ${data[RD_DIR]}".execute().waitFor()
//    }
    
//    Collections.shuffle(failsToDo, new Random())
//    int threads = 5
//    Utils.runCommands(failsToDo.collect{it[RD_COMM]}, threads, failsToDo.collect{it[RD_OUTFILE]}, false)

    println '=============================='
    
    // print commands to stdout
    failsToDo.each{
        println "${it[RD_COMM]} > ${it[RD_OUTFILE]}"
    }
    
    println '=============================='
    
}//end action ACT_SPLIT_FAIL

println "End time: ${new Date()}\n"

return 0