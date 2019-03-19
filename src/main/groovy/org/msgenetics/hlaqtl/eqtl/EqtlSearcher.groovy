/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.msgenetics.hlaqtl.eqtl.LDResult
import java.util.concurrent.locks.ReentrantLock
import java.util.logging.Level
import java.util.logging.Logger
import org.ngsutils.variation.SNPData
import org.msgenetics.hlaqtl.CompletionReporter
import org.ngsutils.Utils

/**
 *
 * @author victor
 */
class EqtlSearcher {
	
    static final String CORREL_FILE = Main.CORREL_FILE+'.filter'
    static final String GROUP_LD_FILE = 'group_ld_result.txt'
    static final String GROUP_CONTROL_LD_FILE = 'group_cntrl_ld_result.txt'
    
    //correlation.txt file fields
    static final int F_SNP = 0
    static final int F_ISO = 1
    static final int F_CORR = 2
    static final int F_PVAL = 5
    static final String SNP_TPED_PREF = SNPManager._1000G_PED
    
    //isoform filter file fields
    static final int ISOF_ID = 0
    static final int ISOF_CHR = 8
    static final int ISOF_START = 9
    static final int ISOF_END = 10
    static final int ISOF_TYPE = 7
    static final int ISOF_ENSGENE = 5
    static final int ISOF_GENE = 6
    
    static final regexTped = /\w+_chr\w+\.tped/
    
    //output files
    static final String EQTL_TRANS_LD = 'eqtls_trans_ld.txt'
    
    double pValThr = 1e-5 // p-value threshold for best-snp search
    double ldR2Thr = 0.4 // LD R-sq threshold used to accept eQTLs
    private File tfamFile = null
    def searchList = [] as TreeSet // snps or isoforms to search in eqtls
    def resultMap = [:] as TreeMap // map with key=(snp or iso) and value=list of associated iso/snps
    def bestEqtls = [:] as TreeMap // best snp of each isoform (key=isoform, value=snp map info)
    def snpData = [:] as TreeMap // map with key=snp and value=SNPData object
    def chrSnps = [:] as TreeMap // map with key=chr and value=list of SNPData objects
    
    
    /**
     *
     */
    public EqtlSearcher() {
        
    }
    
    public EqtlSearcher(Collection list) {
        searchList.addAll(list)
    }
    
    
    /**
     *
     *  @return a Map with {key=snp{ and {value=list of isoforms where snp is true eqtl}
     */
    public static calcTrueEqtls(String listFile, String eqtlsDir, String workDir){
        def snplist = []
        new File(listFile).eachLine{ snplist << it }
        
        EqtlSearcher instance = new EqtlSearcher(snplist)
        
        instance.search(eqtlsDir, F_SNP)
        
        // get isoforms related to search snps
        def isoSet = instance.resultMap.collect{ it.value }.flatten() as TreeSet
        
        instance.searchBestEqtls(eqtlsDir, isoSet)
        
        instance.generateTpedFiles(workDir)
        
        def trueEqtls = [:]
        
        instance.resultMap.each{ snp, isolist->
            trueEqtls[snp] = []
            def isolist2 = isolist.findAll{ instance.bestEqtls.containsKey(it) }
            
            isolist2.each{ iso->
                def bestSnp = instance.getBestEqtl(iso)
                
                try{
                    if( bestSnp==snp || instance.testLD(snp,bestSnp) ){
                        trueEqtls[snp] << iso
                    }
                } 
                catch(ex){
                    Logger.getLogger(EqtlSearcher.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        return trueEqtls
    }
    
    /**
     * parse groups file (format: group snp_id)
     * 
     */ 
    static protected parseSnpGroups(file, Map groups, snplist) {
        def reader = new File(file).newReader()
        reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            if( toks[1].startsWith('rs') ) {
                if( !groups.containsKey(toks[0]) ){
                    groups[toks[0]] = [] as TreeSet
                }

                groups[toks[0]] << toks[1]
                snplist << toks[1]
            }
        }
        
        reader.close()
    }
    
    /**
     * performs LD calculation among SNPs in the same group
     */
    public static performGroupLD(args){
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String groupsFile = Main.readOption(args, Main.OPT_SNPS)
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        
        def dirs = [eqtlCisDir,outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        eqtlCisDir = dirs[0]
        outdir = dirs[1]
        
        println "output: ${outdir}"
        println "Groups file: ${groupsFile}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        
        def groups = [:] as TreeMap
        def snplist = []
        
        //parse groups file (format: group snp_id)
        parseSnpGroups(groupsFile, groups, snplist)
        
        //search SNPs in eqtl cis result
        EqtlSearcher instance = new EqtlSearcher(snplist)
        instance.search(eqtlCisDir, F_SNP)
        instance.generateTpedFiles(outdir)
        
        //perform LD calculations and write output
        def writer = new PrintWriter(outdir+GROUP_LD_FILE)
        writer.println("group\tsnp1\tsnp2\trSquare")
        
        groups.each{ name, set->
            def list = set as List
            for(i in 0..list.size()-2){
                for(j in i+1..list.size()-1){
                    Double rSq = null
                    try{
                        rSq = LDCalculator.calcLD(instance.getSnpData(list[i]), instance.getSnpData(list[j]))
                    } catch(e){ }
                    writer.println("${name}\t${list[i]}\t${list[j]}\t${(rSq==null) ? 'Error' : rSq}")
                }
            }
        }
        
        writer.close()
    }
    
    /**
     * performs LD calculation among SNPs in the control group and SNPs in each others
     */
    public static performControlGroupLD(args) {
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String groupsFile = Main.readOption(args, Main.OPT_SNPS)
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        
        def dirs = [eqtlCisDir,outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        eqtlCisDir = dirs[0]
        outdir = dirs[1]
        
        println "output: ${outdir}"
        println "Groups file: ${groupsFile}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        
        def groups = [:] as TreeMap
        def snplist = [] as TreeSet
        
        //parse groups file (format: group snpId chr position)
        parseSnpGroups(groupsFile, groups, snplist)
        
        //search SNPs in eqtl cis result
        EqtlSearcher instance = new EqtlSearcher(snplist)
        
        if( !instance.existsTped(outdir) ){
            instance.searchSnps(eqtlCisDir, snplist)
            instance.generateTpedFiles(outdir)
        }
        else{
            println "tped files already exists. loading from ${outdir}"
            instance.loadTpedFiles(outdir)
        }
        
        //perform LD calculations
        def noncontrol = groups.keySet().findAll{ it!='control' }
        def control = groups['control']
        def writerLock = new ReentrantLock() 
        
        
        // build map for LD results
        def buildKey = { snp1, snp2-> Utils.hash(snp1+':'+snp2) }
        def ldRsqResults = [:] as TreeMap
        control.each{ ldRsqResults[it] = [:] as TreeMap }
        
        int nonCtrlSnps = noncontrol.sum{groups[it].size()}
        
        def reporter = new CompletionReporter( title:'LD calc', 
            totalJobs:(nonCtrlSnps * control.size()) )
        
        //closure for LD calculation (threaded)
        def lDThreadCalc = { thread, nthreads ->
            def thrLdRsqResults = [:] as TreeMap
            
            control.eachWithIndex{ snp1, i->
                if( (i%nthreads)==thread ) {
                    def resMap = thrLdRsqResults[snp1]
                    if( resMap==null ){
                        resMap = [:] as TreeMap
                        thrLdRsqResults[snp1] = resMap
                    }
                    
                    def data1 = instance.getSnpData(snp1)
                    
                    noncontrol.each{ name->
                        def grp = groups[name]
                        
                        grp.each{ snp2->
                            if( !resMap.containsKey(snp2) ){
                                Double rSq = null
                                def data2 = instance.getSnpData(snp2)

                                try{
                                    if( data1 && data2 && data1.chrNum==data2.chrNum ) {
                                        rSq = LDCalculator.calcLD(data1, data2)
                                    }
                                } catch(ex) { 
                                    System.err.println("Error in LD calc. ${snp1}-${snp2}: ${ex.message}")
                                }

                                // store result
                                if( rSq!=null ){ resMap[snp2] = rSq }
                            }
                        }
                    }
                    
                    reporter.updateSafe( nonCtrlSnps )
                }
            }
            
            writerLock.lock()
            // add thread results to instance results
            thrLdRsqResults.each{ snp, map->
                ldRsqResults[snp] += map
            }
            writerLock.unlock()
        }//end closure lDThreadCalc
        
        //run lDThreadCalc in threads
        int nthreads = 6
        Utils.runClosures(
            [{lDThreadCalc(0, nthreads)}, {lDThreadCalc(1, nthreads)},
            {lDThreadCalc(2, nthreads)}, {lDThreadCalc(3, nthreads)},
            {lDThreadCalc(4, nthreads)}, {lDThreadCalc(5, nthreads)}/*,
            {lDThreadCalc(6, nthreads)}, {lDThreadCalc(7, nthreads)}*/], 
            nthreads )
        
        // write output
        def writer = new PrintWriter(outdir+GROUP_CONTROL_LD_FILE)
        writer.println("group\tsnp_control\tlocus1\tsnp2\tlocus2\trSquare")
        
        control.each{ snp1->
            def data1 = instance.getSnpData(snp1)
            
            noncontrol.each{ name->
                def grp = groups[name]
                grp.each{ snp2->
                    def data2 = instance.getSnpData(snp2)
                    Double rSq = ldRsqResults[snp1][snp2]
                    
                    def result = rSq
                    boolean skip = false
                    
                    if( data1==null || data2==null ){
                        result = 'missing_snp'
                    }
                    else if( data1 && data2 && data1.chrNum!=data2.chrNum ){
                        skip = true
                    }
                    
                    if( !skip ) {
                    writer.println("${name}\t${snp1}\t${data1?.locus}\t${snp2}\t${data2?.locus}\t${result}")
                    }
                }
            }
        }
        
        writer.close()
    }
    
    /**
     * generate tped files for given snps
     * 
     */ 
    static generateGenotypes(args) {
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String snpsFile = Main.readOption(args, Main.OPT_SNPS)
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        
        def dirs = [eqtlCisDir,outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        eqtlCisDir = dirs[0]
        outdir = dirs[1]
        
        println "output: ${outdir}"
        println "snps file: ${snpsFile}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        
        // get snps from file
        def snplist = []
        new File(snpsFile).eachLine{ line-> snplist << line }
        
        EqtlSearcher instance = new EqtlSearcher(snplist)
        instance.search(eqtlCisDir, F_SNP)
        instance.generateTpedFiles(outdir)
    }
    
    /**
     * performs actual eqtl-trans calculation. Calculates LD between a eqtl-trans SNP
     * and the best cis SNP for the eqtl-trans isoform.
     * 
     */
    static performActualEqtlTransLD(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(Main.OPT_INPUT_DIR) }
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String isoFile = Main.readOption(args, Main.OPT_ISOFORMS)
        
        def dirs = [outdir, inputdir, eqtlCisDir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        inputdir = dirs[1]
        eqtlCisDir = dirs[2]
        
        println "output: ${outdir}"
        println "input eqtls-trans dir: ${inputdir}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        println "isoforms filtered file: ${isoFile}"
        
        def isoformsData = [:] as TreeMap
        
        // 1 - parse isoforms filtered file
        def reader = new File(isoFile).newReader()
        reader.readLine() //skip header
        
        reader.splitEachLine("\t"){ toks->
            isoformsData[toks[0]] = buildIsoData(toks)
        }
        
        reader.close()
        
        // 2 - get best eqtls cis for selected isoforms
        EqtlSearcher instance = new EqtlSearcher()
        instance.pValThr = 0.05
        instance.searchBestEqtls(eqtlCisDir, isoformsData.keySet())
        
        // 3 - read eqtl-trans and generate Tped files
        EqtlSearcher transSearch = new EqtlSearcher()
        def isoTransCorrMap = transSearch.searchIsoformsExtended(inputdir, isoformsData.keySet())
        
        def snpSet = [] as TreeSet
        isoTransCorrMap.each{ iso, list-> list.each{snpSet << it[F_SNP]} }
        instance.searchSnps(inputdir, snpSet)
        instance.generateTpedFiles(outdir)
        
        // 4 - for each isoform perform LD calculation between best SNP cis and SNPs trans
        def reporter = new CompletionReporter( title:'LD calc', 
            totalJobs:isoTransCorrMap.values().sum{it.size()} )
        
        // log LD results
        def writer = new PrintWriter(outdir+EqtlStatistics.PLINKLD_LOG)
        writer.println("snp1\tlocus1\tsnp2\tlocus2\trSquare")
        
        def snpsNotFound = [] as Set
        def writerLock = new ReentrantLock() 
        
        // build map for LD results
        def ldRsqResults = [:] as TreeMap
        isoTransCorrMap.each{ isoId, list->
            def id = instance.bestEqtls[isoId]?.get(F_SNP)
            
            if( id!=null && !ldRsqResults.containsKey(id) ){
                ldRsqResults[id] = [:] as TreeMap
            }
        }
        
        //closure for LD calculation (threaded)
        def lDThreadCalc = { thread, nthreads ->
            def thrLdRsqResults = [:] as TreeMap
            def thrSnpsNotFound = [] as Set
            
            isoTransCorrMap.eachWithIndex{ isoId, list, i->
                if( (i%nthreads)==thread && (instance.bestEqtls[isoId]?.get(F_SNP)!=null) ){
                    def data = instance.getSnpData(instance.bestEqtls[isoId][F_SNP])
                    
                    def resMap = thrLdRsqResults[data.id]
                    if( resMap==null ){
                        resMap = [:] as TreeMap
                        thrLdRsqResults[data.id] = resMap
                    }
                    
                    StringBuilder logBuffer = new StringBuilder()

                    list.each{ corrObj->
                        try{
                            def data2 = instance.getSnpData(corrObj[F_SNP])
                            if( data2!=null && !resMap.containsKey(data2.id) ){
                                Double rSq = null

                                if( data.chrNum==data2.chrNum) {
                                    rSq = (data.id==data2.id) ? 1.0 : LDCalculator.calcLD(data, data2)
                                }

                                // store result
                                if( rSq!=null ){
                                    resMap[data2.id] = new LDResult(snpTrans:data2, snpBestCis:data, rSq:rSq)
                                }

                                // log LD result
                                logBuffer << "${data.id}\t${data.locus}\t${data2.id}\t${data2.locus}\t${(rSq==null) ? 'NA' : rSq}\n"
                            }
                            else if( data2==null ){
                                thrSnpsNotFound << data2.id
                            }
                        } catch(ex) {
                            System.err.println("Error in LD calc. ${instance.bestEqtls[isoId][F_SNP]}-${corrObj[F_SNP]}: ${ex.message}")
                        }
                    }

                    reporter.updateSafe(list.size())
                    writerLock.lock()
                    writer.print(logBuffer.toString())
                    writerLock.unlock()
                }
            }
            
            writerLock.lock()
            // add thread results to instance results
            thrLdRsqResults.each{ snp1, map-> ldRsqResults[snp1] += map }
            snpsNotFound.addAll(thrSnpsNotFound)
            writerLock.unlock()
        }//closure lDThreadCalc
        
        //run lDThreadCalc in threads
        int nthreads = 8
        Utils.runClosures(
            [{lDThreadCalc(0, nthreads)}, {lDThreadCalc(1, nthreads)},
            {lDThreadCalc(2, nthreads)}, {lDThreadCalc(3, nthreads)},
            {lDThreadCalc(4, nthreads)}, {lDThreadCalc(5, nthreads)},
            {lDThreadCalc(6, nthreads)}, {lDThreadCalc(7, nthreads)}], 
            nthreads )
        
        writer.close()
        
        //print snps not found
        snpsNotFound.each{ System.err.println("SNP ${it} not found.") }
        
        // End - write result
        writer = new PrintWriter(outdir+EQTL_TRANS_LD)
        writer.println("snp\tsnp_locus\tiso\tiso_locus\tbiotype\tcorr\tpval\tsnp_cis\tsnp_cis_locus\tLD_rSq\ttype")//header
        
        isoTransCorrMap.each{ isoId, list->
            def snpCis = instance.bestEqtls[isoId]?.get(F_SNP)
            def data = (snpCis!=null) ? instance.getSnpData(snpCis) : null
            def isoData = isoformsData[isoId]
            def isoLocus = "${isoData[ISOF_CHR]}:${isoData[ISOF_START]}-${isoData[ISOF_END]}"
            
            list.each{ corrObj->
                def data2 = instance.getSnpData(corrObj[F_SNP])
                LDResult res = (data!=null) ? ldRsqResults[data.id][data2.id] : null
                def type = (isEqtlTrans(data2,isoData) ? 'trans' : 'cis')
                
                writer.println(
"${data2.id}\t${data2.locus}\t${isoId}\t${isoLocus}\t${isoData[ISOF_TYPE]}\t${corrObj[F_CORR]}\t${corrObj[F_PVAL]}\t${data?.id}\t${data?.locus}\t${res?.rSq}\t${type}"
                )
            }
        }
        
        writer.close()
    }
    
    
    /**
     * used by performActualEqtlTransLD
     */
    static private buildIsoData(toks) {
        return [ (ISOF_ID):toks[ISOF_ID], (ISOF_CHR):toks[ISOF_CHR], 
            (ISOF_START):toks[ISOF_START] as Integer, (ISOF_END):toks[ISOF_END] as Integer, 
            (ISOF_TYPE):toks[ISOF_TYPE], (ISOF_ENSGENE):toks[ISOF_ENSGENE], 
            (ISOF_GENE):toks[ISOF_GENE] ]
    }
    
    /**
     * used by performActualEqtlTransLD
     */
    static boolean isEqtlTrans(SNPData snp, Map iso) {
        if( iso[ISOF_CHR]!=snp.chrNum ){
            return true
        }    
        double dist = Math.min( Math.abs(snp.position-iso[ISOF_START]), Math.abs(snp.position-iso[ISOF_END]) )
        return (dist > EqtlStatistics.TRANS_EQTL_THR)
    }
    
    /**
     * search eqtl data for snps/isoforms in this.searchList 
     */
    def search(dir, fieldIdx) {
        int resIdx = (fieldIdx==F_SNP) ? F_ISO : F_SNP
        
        new File(dir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE)

            if( file.exists() ){
                def reader = file.newReader()
                reader.readLine()
                boolean found = false
                def snpCurrList = [] as TreeSet

                reader.splitEachLine("\\s"){ toks->
                    def field = toks[fieldIdx]
                    
                    if( field in searchList ){
                        if( resultMap[field]==null )
                            resultMap[field] = []
                            
                        found = true
                        resultMap[field] << toks[resIdx]
                        
                        if(fieldIdx==F_ISO){
                            snpCurrList << toks[resIdx]
                        }
                    }
                }

                reader.close()
                
                if(found) {
                    // get snp info from .tped files
                    getTpedSnpInfo(fDir, ((fieldIdx==F_SNP) ? searchList : snpCurrList))
                }
            }
        }
    }
    
    /**
     *
     */
    def searchPVal(dir, fieldIdx, double pVal) {
        int resIdx = (fieldIdx==F_SNP) ? F_ISO : F_SNP
        
        new File(dir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE)

            if( file.exists() ){
                def reader = file.newReader()
                reader.readLine()
                
                reader.splitEachLine("\\s"){ toks->
                    def field = toks[fieldIdx]
                    double value = toks[F_PVAL] as Double
        
                    if( value<pVal ) {
                        if( resultMap[field]==null )
                            resultMap[field] = []
                            
                        resultMap[field] << toks[resIdx]
                    }
                }
                
                reader.close()
            }
        }
    }
    
    /**
     *
     */
    private void getTpedSnpInfo(fDir, snpSet) {
        // get snp info from .tped files
        def tpedFiles = (fDir.listFiles() as List).findAll{ 
            it.name.startsWith(SNP_TPED_PREF+'_chr') && it.name.endsWith('.tped') }
        
        if( !this.tfamFile ) {
            this.tfamFile = (fDir.listFiles() as List).find{ 
                it.name.startsWith(SNP_TPED_PREF+'_chr') && it.name.endsWith('.tfam') }
        }

        if( !tpedFiles ) {
            tpedFiles = (fDir.listFiles() as List).findAll{ 
                it.name==(SNP_TPED_PREF+'.tped') }
            
            if( !this.tfamFile ) {
                this.tfamFile = (fDir.listFiles() as List).find{ 
                    it.name==(SNP_TPED_PREF+'.tfam') }
            }
        }
                    
        tpedFiles?.each{ tpedFile->
            def reader = tpedFile.newReader()

            reader.eachLine{ line->
                def toks = line.split("\\s",5)
                def id = toks[1]
                if( (id in snpSet) && !snpData.containsKey(id) ){
                    // add snp info to snpData map
                    def data = new SNPData(id:id, chr:toks[0], 
                            position:(toks[3] as Integer), reference:tpedFile,
                            genotypes:toks[4])
                        
                    snpData[id] = data
                        
                    if( chrSnps[data.chrNum]==null ) {
                        chrSnps[data.chrNum] = []
                    }
                    
                    chrSnps[data.chrNum] << data
                }
            }

            reader.close()
        }
    }
    
    /**
     * finds the best eqtl of each isoform in isoSet
     */
    def searchBestEqtls(dir, isoSet) {
        new File(dir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE)

            if( file.exists() ){
                searchBestEqtlsInternal(file, isoSet, fDir)
            }
        }
    }
    
    /**
     * collect snp data for snps in snpList
     */
    def searchSnps(dir, snpList){
        new File(dir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE)

            if( file.exists() ){
                getTpedSnpInfo(fDir, snpList)
            }
        }
    }
    
    /**
     * collect isoform extended correlation data for isoforms in isoSet
     * 
     * @return a map with key=isoformId and value=list of eqtl data
     */
    def searchIsoformsExtended(dir, isoSet){
        def isoCorrMap = [:] as TreeMap
        
        new File(dir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE)

            if( file.exists() ){
                def snpSet = parseCorrFile(file, isoSet, isoCorrMap)
                getTpedSnpInfo(fDir, snpSet)
            }
        }
        
        return isoCorrMap
    }
    
    /**
     *
     * @return a Set of present eQTLs snps
     */
    private def parseCorrFile(File file, isoSet, isoCorrMap, double pPval = 0.05) {
        def reader = file.newReader()
        reader.readLine()//skip header
        
        def snpSet = [] as TreeSet
        
        reader.splitEachLine("\\s"){ toks->
            double pval = toks[F_PVAL] as Double

            if( pval<pPval && (toks[F_ISO] in isoSet) ){
                //store snp-iso correlation value
                if( isoCorrMap[toks[F_ISO]]==null )
                    isoCorrMap[toks[F_ISO]] = []

                isoCorrMap[toks[F_ISO]] << [ (F_ISO):toks[F_ISO], (F_SNP):toks[F_SNP],
                    (F_CORR):(toks[F_CORR] as Double), (F_PVAL):pval ]
                
                snpSet << toks[F_SNP]
            }
        }

        reader.close()
        
        return snpSet
    }
    
    /**
     *
     */
    private def searchBestEqtlsInternal(File file, isoSet, fileDir) {
        def isoCorrMap = [:]
        parseCorrFile(file, isoSet, isoCorrMap, pValThr)

        //get the best snp of each isoform
        def bestSnpSet = [] as TreeSet
        
        isoCorrMap.each{ iso, list->
            def bestSnp = list.max{ Math.abs(it[F_CORR]) }
            bestSnpSet << bestSnp[F_SNP]
            
            def current = bestEqtls[iso]
            
            if( current==null || bestSnp[F_CORR]>current[F_CORR] )
                bestEqtls[iso] = bestSnp
        }
        
        // get tped SNP info for the best SNPs
        getTpedSnpInfo(fileDir, bestSnpSet)
    }
    
    
    /**
     *
     */
    private def generateTpedFiles( String workDir ) {
        
        chrSnps.each{ chr, list->
            //generate tped file for each chromosome
            Collections.sort(list, SNPData.comparator)
            
            def writer = new PrintWriter( compTpedName(workDir, chr) )
            list.each{ writer.println(it.tpedLine) }
            writer.close()
        }
    }
    
    /**
     * composes the generated tped file name for dir and chr number
     */
    private static String compTpedName(dir, chrNum) {
        return "${dir}${SNP_TPED_PREF}_chr${chrNum}.tped"
    }
    
    /**
     * 
     */ 
    protected boolean existsTped(dir) {
        def files = new File(dir).list({d, f-> f ==~ regexTped } as FilenameFilter).toList()
        return !files.isEmpty()
    }
    
    /**
     * load snpData from tped files
     */
    protected loadTpedFiles(dir) {
        def files = new File(dir).list({d, f-> f ==~ regexTped } as FilenameFilter).toList()
        
        files.each{ tpedFile->
            SNPData.createFromTped(dir+tpedFile, searchList, snpData)
        }
    }
    
    /**
     * run LD test for a pair of SNPs
     * 
     * @return true if LD R-sq > ldR2Thr
     */
    private boolean testLD(String snp1, String snp2) {
        SNPData data1 = snpData[snp1]
        SNPData data2 = snpData[snp2]
        
        if( data1.chrNum!=data2.chrNum ){
            return false
        }
        
        Double rSq = LDCalculator.calcLD(data1, data2)
        
        if( rSq==null ) {
            throw new RuntimeException("Error computing R-sq")
        }
        
        return (rSq>this.ldR2Thr)
    }
    
    /**
     *
     */
    public String getBestEqtl(String isoform) {
        def map = bestEqtls[isoform]
        return ((map) ? map[F_SNP] : null)
    }
    
    /**
     *
     */
    public SNPData getSnpData(String snp) {
       return snpData[snp]
    }
    
    /**
     * calculates SNPs data frequencies (MAF, minor/major allele)
     * 
     */
    public def calSnpFrequencies() {
        snpData.each{ id, data-> data.freqA1 }
    }
    
    /**
     * splits the given string and returns the tokens which are not spaces/blank
     * characters (' ', '\t', ...)
     */
    private static getClearTokens(String string){
        return (string.split("\\s") as List).findAll{!it.isEmpty()}
    }
    
}

