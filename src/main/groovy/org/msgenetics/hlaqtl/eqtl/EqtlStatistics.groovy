/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.ngsutils.EnsemblUtils
import org.msgenetics.hlaqtl.SimpleBLAST
import org.msgenetics.hlaqtl.BLASTResult
import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import org.ngsutils.semantic.query.GeneQueryUtils
import org.ngsutils.semantic.LinkedLifeDataFactory
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import org.jfree.data.statistics.SimpleHistogramDataset
import org.jfree.data.statistics.SimpleHistogramBin
import org.jfree.chart.JFreeChart
import org.jfree.chart.plot.XYPlot
import org.jfree.chart.plot.PlotOrientation
import org.jfree.chart.axis.NumberAxis
import org.jfree.chart.renderer.xy.XYBarRenderer
import org.jfree.chart.ChartFactory
import org.jfree.chart.ChartUtilities
import java.util.logging.Level
import java.util.logging.Logger
import org.msgenetics.hlaqtl.CompletionReporter

/**
 *
 */
class LDResult {
    SNPData snpBestCis
    SNPData snpTrans
    SNPData snpFeat = null
    Double rSq
}//

/**
 *
 */
class SNPFeatures {
    static final MISSING_CODES = ['','NA','NR','-','null']
    
    def snps = [:] as TreeMap //key=snp, value=set of features
    def features = [:] as TreeMap //key=feature, value=set of snps
    
    /**
     * file format: <chr>\t<snps>\t<features/diseases>
     */
    static SNPFeatures loadFeatSnps(String file) {
        def reader = new File(file).newReader()
        reader.readLine() //skip header
        
        SNPFeatures featObj = new SNPFeatures()
        
        reader.splitEachLine("\t"){ toks->
            def chr = toks[0]
            def snps = toks[1]
            def feats = toks[2]
            featObj.addFeatures(chr, snps, feats)
        }
        
        reader.close()
        
        return featObj
    }
    
    /**
     *
     */
    protected def splitField(str){
        if( str in MISSING_CODES ){
            return null
        }
        return (str.split(",") as List).collect{it.trim()}
    }
    
    /**
     *
     */
    protected def addToMap(id, Map map, list){
        if( !map.containsKey(id) ){
            map[id] = [] as Set
        }
        map[id].addAll(list)
    }
    
    /**
     *
     */
    protected def addFeatures(String chr, String snpsF, String featsF){
        def listSnps = splitField(snpsF)
        def listFeats = splitField(featsF)
        
        if( listSnps && listFeats ){
            listSnps.each{snp-> addToMap(snp,snps,listFeats) }
            listFeats.each{feat-> addToMap(feat,features,listSnps) }
        }
    }
    
    /**
     *
     */
    Set getSnpSet(){
        return snps.keySet()
    }
    
    /**
     *
     */
    Set getSnpFeatures(id){
        return snps[id]
    }
    
    /**
     *
     */
    Set getFeatureSet(){
        return features.keySet()
    }
    
}//end of SNPFeatures


/**
 *
 * @author victor
 */
class EqtlStatistics {
	
    static final String CORREL_FILE_UNFILT = Main.CORREL_FILE
    static final String SNP_TPED_PREF = SNPManager._1000G_PED
    static final String SEQS_DIR = 'seqs'
    static final String TPED_DIR = 'tped'
    static final String TAX_ID = '9606'
    static final String CORR_FILE_COLS = "snp\tisoform\tcorrelation\tpvalue\tlog_pvalue\tcorr_pvalue\tlog_corr_pvalue\tnum_obs\ttype"
    static final String CORR_FILE_EXTCOLS = "\tsnp_locus\tisoform_locus\tiso_eqtl_cis\tblast_ident\tblast_score"
    static final String CORR_FILE_LDCOLS = "\tbest_snp_iso_cis\tbest_snp_cis_locus\tld_rsquare"
    
    static final String BLAST_LOG_HEAD = "snp\tisoform_cis\tisoform_trans\t#alignment\tscore_bits\tscore"+
        "\texpect\tidentities\tgaps\tstrand\tquery\tlen\tsubject\tlen"
    
    // log files to generate
    static final String BLAST_LOG = 'blastn.log'
    static final String PLINKLD_LOG = 'plink_ld.log'
    static final String ISO_EQTL_CIS_TRANS = 'iso_eqtls_cis_trans.txt'
    static final String TRAIT_EQTL_ASSOC = 'trait_eqtl_assoc.txt'
    
    //correlation.txt file fields
    static final int F_SNP = 0
    static final int F_ISO = 1
    static final int F_CORR = 2
    static final int F_PVAL = 3
    static final int F_LOG_PVAL = 4
    static final int F_FDR_PVAL = 5
    static final int F_LOG_FDR_PVAL = 6
    static final int F_NOBS = 7
    
    //another constants
    static final double BIN_STEP = 1.0/0.25
    static final double TRANS_EQTL_THR = 1e6 //inter-chromosome trans eqtl distance threshold
    static final int CHART_WIDTH = 1280
    static final int CHART_HEIGHT = 1024
    
    Map bins = [:] as TreeMap //key=bin index, value=count of eqtls with p-value in bin
    Map binsTrans = [:] as TreeMap
    Map isoforms = [:] as TreeMap //key=isoform_id, value IsoformExpData object
    Map snpData = [:] as TreeMap //key=snp_id, value SNPData object
    Map transEqtls = [:] as TreeMap //key=snp_id, value=list of isoforms_ids
    Map blastResults = [:] as TreeMap //key="iso_cis_id:iso_trans_id", value=BLASTResult object
    Map ldRsqResults = [:] as TreeMap //key=snp_id, value=map{key=isoform_id, value=LDResult Object}
    
    def graph //semantic data graph
    GeneQueryUtils geneQuery
    String eqtlDir
    SummaryStatistics summary
    SummaryStatistics summaryTrans
   
    //filtering test
    public enum FilterTest {CORR_PVAL, PERCENTILE}
    FilterTest filter = FilterTest.PERCENTILE
    
    //closure 
    static blastResKey = {isoCis, isoTrans-> "${isoCis}:${isoTrans}"}
    
    /**
     *
     */
    public EqtlStatistics(String p_eqtlDir) {
        eqtlDir = p_eqtlDir
        parseSnpAndIsoforms()
    }
    
    
    /**
     * perform eqtl-trans filtering
     * 
     */
    static perform(args) {
        //read input and output dir
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(Main.OPT_INPUT_DIR) }
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String percStr = Main.readOption(args, Main.OPT_VALUE)
        String blast = Main.readOption(args, Main.OPT_BLAST)
        String cisEqtls = Main.readOption(args, Main.OPT_SNPS)
        
        def dirs = [outdir, inputdir, eqtlCisDir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        inputdir = dirs[1]
        eqtlCisDir = dirs[2]
        
        println "output: ${outdir}"
        println "input eqtls dir: ${inputdir}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        println "cis eqtls file: ${cisEqtls}"
        println "BLAST+ path: ${blast}"
        
        ///double percentile = (!percStr) ? 95.0 : percStr as Double
        double pvalCutoff = (!percStr) ? 0.05 : percStr as Double
        boolean transOnly = true
        
        EqtlStatistics instance = new EqtlStatistics(inputdir)
        instance.calcSummary()
        instance.filter = FilterTest.CORR_PVAL
        instance.loadSemantic(outdir)
        instance.generateStatistics(outdir)
        instance.generateCorrFile(outdir, pvalCutoff, transOnly)
        
        
        instance.parseCisEqtls(cisEqtls, 1, 2)//snpIdx=1, corrIdx=2
        instance.searchBestSnpsCis(outdir, eqtlCisDir)
        instance.downloadSeqs(outdir) // download ensembl sequences
        instance.compareSeqs(outdir, blast) // perform blast among cis and trans eqtl isoforms
        instance.generateCorrFileExtended(outdir)
    }
    
    
    /**
     * perform search of isoforms that are eqtls cis and trans
     * 
     */
    static performSearchCisTrans(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(Main.OPT_INPUT_DIR) }
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        
        def dirs = [outdir, inputdir, eqtlCisDir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        inputdir = dirs[1]
        eqtlCisDir = dirs[2]
        
        double pValCis = 1e-6
        double pValTrans = 0.05
        
        println "output: ${outdir}"
        println "input eqtls-trans dir: ${inputdir}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        
        
        // search eqtl-cis
        def searcherCis = new EqtlSearcher()
        searcherCis.searchPVal(eqtlCisDir, EqtlSearcher.F_ISO, pValCis)
        def isoformsCis = searcherCis.resultMap.keySet()
        
        println "${isoformsCis.size()} eqtls-cis below pvalue ${pValCis}. ${new Date()}"
        
        EqtlStatistics instance = new EqtlStatistics(inputdir)
        instance.collectTransEqtls(pValTrans)
        
        // get eqtl-trans that are also eqtl-cis
        def isoformsCisTrans = [:] as TreeMap
        
        instance.transEqtls.each{ snpId, list->
            list.each{ isoId->
                if( isoId in isoformsCis) {
                    if( !isoformsCisTrans.containsKey(isoId) ){
                        isoformsCisTrans[isoId] = [] as Set
                    }
                    
                    isoformsCisTrans[isoId] << snpId
                }
            }
        }
        
        println "${isoformsCisTrans.size()} eqtls trans/cis below pvalue ${pValTrans}. ${new Date()}"
        
        //write output
        def outFile = ISO_EQTL_CIS_TRANS
        def writer = new PrintWriter(outdir+outFile)
        writer.println("isoform\tlocus\tsnps_cis\tsnps_trans") 
        
        isoformsCisTrans.each{ isoId, snpsTrans->
            def iso = instance.isoforms[isoId]
            String locus = "${iso.chr}:${iso.start}-${iso.end}"
        
            writer.println("${isoId}\t${locus}\t${searcherCis.resultMap[isoId]}\t${snpsTrans}") 
        }
        
        writer.close()
    }
    
    /**
     * Builds a table containing associations between eqtl snps and featured snps (disease).
     * The association is calculated using LD rSquare.
     * 
     */
    static performAssocEqtlSnps(args){
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(Main.OPT_INPUT_DIR) }
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String cisEqtls = Main.readOption(args, Main.OPT_SNPS)
        String featSnps = Main.readOption(args, Main.OPT_VALUE)
        
        def dirs = [outdir, inputdir, eqtlCisDir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        inputdir = dirs[1]
        eqtlCisDir = dirs[2]
        
        println "output: ${outdir}"
        println "input eqtls-trans dir: ${inputdir}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        println "cis eqtls file: ${cisEqtls}"
        println "featured SNPs file: ${featSnps}"
        
        def snplist = []
        ///double pValTrans = 0.05
        
        // load featured SNPs file (SNPs, chr and associated terms/features/diseases)
        SNPFeatures featObj = SNPFeatures.loadFeatSnps(featSnps)
        snplist.addAll( featObj.getSnpSet() )
        
        println "${featObj.snps.size()} trait/disease snps loaded"
        println "${featObj.features.size()} traits/diseases loaded"
        
        // load best eqtls and associated trans eqtls
        EqtlStatistics instance = new EqtlStatistics(inputdir)
        ///instance.collectTransEqtls(pValTrans)
        instance.parseCisEqtls(cisEqtls, 1, 2)//snpIdx=1, corrIdx=2
        println "${instance.snpData.size()} eqtl snps loaded"
        snplist.addAll( instance.snpData.keySet() )
        
        //collect snp data
        EqtlSearcher searcher = new EqtlSearcher()
        searcher.searchSnps(eqtlCisDir, snplist)
        searcher.generateTpedFiles(outdir)
        
        // calculate LD between eqtl snps and features snp
        def reporter = new CompletionReporter( title:'LD calc', 
            totalJobs:instance.snpData.size()*featObj.snps.size() )
        
        // log LD results
        def writer = new PrintWriter(outdir+PLINKLD_LOG)
        writer.println("snp1\tlocus1\tsnp2\tlocus2\trSquare")
        
        def snpsNotFound = [] as Set
        def writerLock = new Object()
        
        def lDThreadCalc = { thread, total ->
            def thrLdRsqResults = [:] as TreeMap
            def thrSnpsNotFound = [] as Set
            
            instance.snpData.eachWithIndex{ snpId, data, i->
                if( (i%total) == thread ){
                    thrLdRsqResults[snpId] = [:]
                    StringBuilder logBuffer = new StringBuilder()

                    featObj.getSnpSet().each{ snpFeatId->
                        try{
                            def data2 = searcher.getSnpData(snpFeatId)
                            if( data2!=null ){
                                Double rSq = null

                                if( data.chrNum==data2.chrNum) {
                                    rSq = (snpId==snpFeatId) ? 1.0 : LDCalculator.calcLD(data, data2)
                                }

                                // store result
                                if( rSq!=null ){
                                    thrLdRsqResults[snpId][snpFeatId] = new LDResult(snpTrans:data, snpFeat:data2, rSq:rSq)
                                }

                                // log LD result
                                logBuffer << "${snpId}\t${data.locus}\t${snpFeatId}\t${data2.locus}\t${(rSq==null) ? 'NA' : rSq}\n"
                            }
                            else{
                                thrSnpsNotFound << snpFeatId
                            }
                        } catch(ex) {
                            Logger.getLogger(EqtlStatistics.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }

                    reporter.updateSafe(featObj.snps.size())
                    synchronized(writerLock) { writer.print(logBuffer.toString()) }
                }
            }
            
            synchronized(writerLock) {
                // add thread results to instance results
                instance.ldRsqResults.putAll(thrLdRsqResults)
                snpsNotFound.addAll(thrSnpsNotFound)
            }
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
        snpsNotFound.each{ System.err.println("trait/disease SNP ${it} not found.") }
        
        // write output table to file
        def outFile = TRAIT_EQTL_ASSOC
        writer = new PrintWriter(outdir+outFile)
        writer.println("snp\tiso_cis\tsnp_trait\ttrait/disease\tLD_rSq") 
        
        instance.ldRsqResults.each{ snpId, map->
            map.each{ snpFeatId, ldResult->
                def traits = featObj.getSnpFeatures(snpFeatId).sum{it+','}
                writer.println("${snpId}\t${ldResult.snpTrans.eqtl}\t${snpFeatId}\t${traits}\t${ldResult.rSq}") 
            }
        }
        
        writer.close()
    }
    
    
    
    /**
     * prints out statistics and generates histogram .png files
     */
    public void generateStatistics(workDir) {
        println "Generating histograms.\n"
        generateHistogram(bins, 'cis/trans p-values frequencies', workDir+'hist_all.png')
        generateHistogram(binsTrans, 'trans p-values frequencies', workDir+'hist_trans.png')
        
        println "Summary statistics for all p-values:"
        println summary.toString()
        
        println "Summary statistics for trans p-values:"
        println summaryTrans.toString()
    }
    
    /**
     * add to trans eqtls map
     */
    private def addToTransEqtls(toks){
        if( transEqtls[toks[F_SNP]]==null ){
            transEqtls[toks[F_SNP]] = []
        }

        // add isoform to eqtl-trans map (key=snp_id, value=list of isoforms_id)
        transEqtls[toks[F_SNP]] << toks[F_ISO]
    }
    
    /**
     * parses a previously generated trans-eqtl filter file; used by 
     * generateCorrFile.
     */
    private void parseTransFilterFile(String file) {
        def reader = new File(file).newReader()
        reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            addToTransEqtls(toks)
        }
        
        reader.close()
    }
    
    /**
     * 
     * @param workDir : where to write result file
     * @param cutOff: percentile scaled from 0 - 100, or corrected p-value
     */
    public void generateCorrFile(workDir, double cutOff, boolean transOnly) {
        
        def suff = (transOnly) ? '.trans.filter' : '.all.filter'
        
        def fileName = workDir+CORREL_FILE_UNFILT+suff
        def threshold = null
        
        if( filter==FilterTest.PERCENTILE ) {
            println "Generating ${fileName} for ${cutOff} percentile"

            def binMap = (transOnly) ? binsTrans : bins

            if( !binMap) {
                println "Cannot generate ${fileName}, no data"
                return
            }
        
            threshold = calcPercentile(binMap, cutOff)
        }
        else {
            println "Generating ${fileName} for ${cutOff} corrected p-value"
            threshold = cutOff
        }
        
        //checks whether the file already exists
        if( transOnly && (new File(fileName).exists()) ) {
            println "Filtered eqtl-trans file ${fileName} already exists. parsing it."
            parseTransFilterFile(fileName)
            return
        }
        
        def writer = new PrintWriter(fileName)
        writer.println( CORR_FILE_COLS )
        long count = 0
        
        new File(eqtlDir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE_UNFILT)

            if( file.exists() ) {
                def reader = new BufferedReader( new FileReader(file) )
                reader.readLine() //skip header

                reader.eachLine{ line->
                    def toks = line.split("\\s")
                    
                    try {
                        if( filteringTest(toks, threshold) ) {
                            boolean isTrans = isEqtlTrans(toks[F_SNP],toks[F_ISO])
                            
                            if( isTrans ){
                                addToTransEqtls(toks)
                            }

                            if( !transOnly || isTrans ) {
                                writer.println(line+((isTrans) ? "\ttrans" : "\tcis"))
                                count++
                            }
                        }
                    } catch(e){}
                }

                reader.close()
            }
        }
        
        writer.close()
        
        println "#${count} eqtls written to ${fileName}"
    }
    
    
    /**
     *
     */
    protected boolean filteringTest(lineToks, threshold) {
        if( filter==FilterTest.PERCENTILE ) {
            double val = lineToks[F_LOG_PVAL] as Double
            int binIdx = getBinIdx(val)
            return (binIdx>=threshold)
        }
        
        //FilterTest.CORR_PVAL
        double val = lineToks[F_FDR_PVAL] as Double
        return (val<=threshold)
    }
    
    
    /**
     * parse cis-eqtls, assigns isoforms to snps
     */
    protected def parseCisEqtls(String file, int snpIdx, int corrIdx) {
        def reader = new File(file).newReader()
        
        reader.splitEachLine("\\s"){ toks->
            //{iso} {snp} {corr. score}...
            def data = snpData[toks[snpIdx]]
            double score = Math.abs(toks[corrIdx] as Double)
            
            if( data && (score > data.eqtlScore) ){
                data.eqtl = toks[0]
                data.eqtlScore = score
            }
        }
        
        reader.close()
    }
    
    /**
     * load entrez gene semantic data
     */
    protected def loadSemantic(outdir) {
        println 'Loading semantic data ...'
        this.graph = LinkedLifeDataFactory.loadRepository(LinkedLifeDataFactory.LIST_EGENE, [TAX_ID], outdir)
        this.geneQuery = new GeneQueryUtils(this.graph)
    }
    
    /**
     * traverses a directory where eqtl results are stored and collects isoforms
     * and snp info
     */
    private def parseSnpAndIsoforms() {
        
        new File(eqtlDir).eachDir{ fDir->
            //parse isoforms
            def file = new File(fDir.absolutePath+'/'+Main.MERG_TRACK_FILE)
            
            if( file.exists() ) {
                IsoformParser.parseIsoforms(file, isoforms)
                
                //parse snps
                getTpedSnpInfo(fDir)
            }
        }
    }
    
    /**
     * get snp info from .tped files in eqtl result dir
     */
    private void getTpedSnpInfo(fDir) {
        def tpedFiles = (fDir.listFiles() as List).findAll{ 
            it.name.startsWith(SNP_TPED_PREF+'_chr') && it.name.endsWith('.tped') }
        

        if( !tpedFiles ) {
            tpedFiles = (fDir.listFiles() as List).findAll{ 
                it.name==(SNP_TPED_PREF+'.tped') }
        }
                    
        tpedFiles?.each{ tpedFile->
            def reader = tpedFile.newReader()

            reader.eachLine{ line->
                def toks = line.split("\\s",5)
                def id = toks[1]
                
                if( !snpData.containsKey(id) ) {
                    // add snp info to snpData map
                    def data = new SNPData(id:id, chr:toks[0], 
                        position:(toks[3] as Integer))

                    snpData[id] = data
                }
            }

            reader.close()
        }
    }
    
    /**
     * collects statistics for all and trans eqtls
     */
    protected def calcSummary() {
        summary = new SummaryStatistics()
        summaryTrans = new SummaryStatistics()
        
         new File(eqtlDir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE_UNFILT)

            if( file.exists() ) {
                def reader = file.newReader()
                reader.readLine() //skip header

                reader.splitEachLine("\\s"){ toks->
                    try {
                        double val = toks[F_LOG_PVAL] as Double
                        summary.addValue(val)
                        addToBin(val, bins)
                        
                        if( isEqtlTrans(toks[F_SNP],toks[F_ISO]) ){
                           summaryTrans.addValue(val)
                           addToBin(val, binsTrans)
                        }
                    } catch(e){}
                }

                reader.close()
            }
        }
    }
    
    /**
     * collect trans eqtls
     */
    def collectTransEqtls(double pVal) {
        
         new File(eqtlDir).eachDir{ fDir->
            def file = new File(fDir.absolutePath+'/'+CORREL_FILE_UNFILT)

            if( file.exists() ) {
                def reader = file.newReader()
                reader.readLine() //skip header

                reader.splitEachLine("\\s"){ toks->
                    try {
                        double val = toks[F_FDR_PVAL] as Double
                        
                        if( val<pVal && isEqtlTrans(toks[F_SNP],toks[F_ISO]) ){
                            addToTransEqtls(toks)
                        }
                    } catch(e){}
                }

                reader.close()
            }
        }
    }
    
    /**
     *
     */
    protected int getBinIdx(double val){
        return Math.floor(val*BIN_STEP)
    }
    
    /**
     *
     */
    protected void addToBin(double val, binMap){
        int binIdx = getBinIdx(val)
        Integer current = binMap[binIdx]
        binMap[binIdx] = (current==null) ? 1 : current+1
    }
    
    /**
     *
     */
    protected boolean isEqtlTrans(String snpId, String isoId) {
        def iso = isoforms[isoId]
        def snp = snpData[snpId]
        
        if( iso.chrNum!=snp.chrNum ){ 
            return true
        }    
        
        return (Math.min(Math.abs(snp.position-iso.start),Math.abs(snp.position-iso.end)) > TRANS_EQTL_THR)
    }
    
    /**
     * generates histogram .png picture using jfreechart
     */
    protected generateHistogram(Map bins, String title, String outName){
        SimpleHistogramDataset simplehistogramdataset = new SimpleHistogramDataset(title)
        
        Integer max = bins.keySet().max()
        double step = 1.0/BIN_STEP
        
        if( max==null ){
            println "Cannot generate histogram [title:${title}] [file:${outName}], no data"
            return
        }
        
        (0..max).each{
            SimpleHistogramBin simplehistogrambin = new SimpleHistogramBin(it*step, (it+1)*step, true, false)
            Integer count = bins[it]
            simplehistogrambin.setItemCount( (count==null) ? 0 : count )
            simplehistogramdataset.addBin(simplehistogrambin)
        }
        
        def xAxisLabel = '-log(p_value)'
        def yAxisLabel = null
        
        JFreeChart jfreechart = ChartFactory.createHistogram(title, xAxisLabel, yAxisLabel, 
            simplehistogramdataset, PlotOrientation.VERTICAL, false, false, false)
        
        XYPlot xyplot = (XYPlot)jfreechart.getPlot();
        xyplot.setForegroundAlpha(0.85F);
        xyplot.setDomainPannable(true);
        xyplot.setRangePannable(true);
        NumberAxis numberaxis = (NumberAxis)xyplot.getRangeAxis();
        numberaxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        XYBarRenderer xybarrenderer = (XYBarRenderer)xyplot.getRenderer();
        xybarrenderer.setDrawBarOutline(false);
        
        ChartUtilities.saveChartAsPNG( new File(outName), jfreechart, CHART_WIDTH, CHART_HEIGHT)
    }
    
    /**
     * 
     * @param perc : percentile (scaled from 0 - 100) 
     * @return : bin index
     */
    protected int calcPercentile(Map bins, double perc) {    
        double  totalPerc = bins.values().sum()*perc*0.01
        double sum = 0.0
        int max = bins.keySet().max()
        
        for(int i=0; i<=max; i++) {
            Integer count = bins[i]
            sum += ((count==null) ? 0 : count)
            
            if( sum > totalPerc ) {
                return i
            }
        }
        
        return max
    }
    
    /**
     *
     */
    protected def downloadSeqs(String outdir) {
        def seqsToDown = [] as Set
        
        // collect isoforms ids to download
        snpData.each{ seqsToDown<<it.value.eqtl }
        transEqtls.each{ seqsToDown.addAll( it.value ) }
        
        File seqDir = new File(outdir+SEQS_DIR)
        Utils.createDir(seqDir.absolutePath)
        
        seqsToDown.each{ isoId->
            if( isoId) {
                def geneUri = geneQuery.getGeneByName(isoId, TAX_ID)
                
                if( geneUri==null ){
                    System.err.println("Gene id not found for transcript: ${isoId}")
                }
                else{
                    String ensGeneId = geneQuery.getEnsemblGene(geneUri)

                    if( ensGeneId!=null ) {
                        EnsemblUtils.downTranscriptSeq(ensGeneId, isoId, seqDir)
                    }
                    else {
                        System.err.println("Ensembl gene id not found for transcript: ${isoId}")
                    }
                }
            }
        }
    }
    
    /**
     * compares cis-eqtl and trans-eqtl isoforms for each trans-eqtl snp using BLAST
     * 
     */
    protected def compareSeqs(String outdir, String blastDir) {
        File seqDir = new File(outdir+SEQS_DIR)
        
        def writer = new PrintWriter(outdir+BLAST_LOG)
        // write .log header
        writer.println(BLAST_LOG_HEAD)
        
        def blast = new SimpleBLAST(blastDir)
        
        Double totalJobs = transEqtls.values().sum{it.size()}
        totalJobs = (totalJobs!=null) ? totalJobs : 0.0
        println "${totalJobs as Long} potential BLAST+ alignments to do. ${new Date()}"
        def reporter = new CompletionReporter(title:'BLAST+ alignment', totalJobs:totalJobs)
        
        transEqtls.each{ snpId, list->
            def data = snpData[snpId]
            if( data.eqtl ){//if snp is associated to cis-eqtl isoform
                String query = seqDir.absolutePath+'/'+EnsemblUtils.transcriptFileName(data.eqtl)

                list.each{ 
                    String subject = seqDir.absolutePath+'/'+EnsemblUtils.transcriptFileName(it)
                    BLASTResult res = null
                    
                    try{
                        res = blast.runBlastn(query, subject)
                    }
                    catch(ex){
                        System.err.println("Error running BLAST for: ${data.eqtl} - ${it}")
                        Logger.getLogger(EqtlStatistics.class.getName()).log(Level.SEVERE, null, ex); 
                    }
                    
                    if(res) {
                       blastResults[blastResKey(data.eqtl,it)] = res 
                    }
                    
                    
                    // log result
                    res?.alignments.eachWithIndex{ ali, i->
                        writer.println(
                        "${snpId}\t${data.eqtl}\t${it}\t${i+1}\t${ali.scoreBits}\t${ali.score}\t${ali.expect}\t${ali.identities}"+
                        "\t${ali.gaps}\t${ali.strand}\t${res?.query}\t${res?.queryLength}\t${res?.subject}\t${res?.subjectLength}"
                        )
                    }
                }
            }
            
            reporter.update( list.size() )
        }
        
        writer.close()
        println "BLAST+ alignments completed. ${new Date()}"
    }
    
    /**
     *
     */
    protected def searchBestSnpsCis(String outdir, String eqtlCisDir) {
        def isoSet = [] as TreeSet
        def snpList = []
        
        transEqtls.each{ snpId, list-> 
            snpList << snpId
            isoSet.addAll(list) 
        }
        
        String workDir = outdir+TPED_DIR+'/'
        Utils.createDir(workDir)
        
        
        def eqtlSearcher = new EqtlSearcher(snpList)
        eqtlSearcher.pValThr = 0.05
        eqtlSearcher.search(eqtlCisDir, EqtlSearcher.F_SNP)
        eqtlSearcher.searchBestEqtls(eqtlCisDir, isoSet)
        eqtlSearcher.generateTpedFiles(workDir)

        // log LD results
        def writer = new PrintWriter(outdir+PLINKLD_LOG)
        writer.println("snp\tlocus\tbest_snp_cis\tlocus\tisoform\trSquare")
        
        //calculate LD
        transEqtls.each{ snp, list-> 
            ldRsqResults[snp] = [:]
            def snpData = eqtlSearcher.getSnpData(snp)
            
            list.each{ iso->
                def bestSnpCis = eqtlSearcher.getBestEqtl(iso)
                def dataBestSnp = (bestSnpCis) ? eqtlSearcher.getSnpData(bestSnpCis) : null
                
                if(dataBestSnp) {
                    try{
                        Double rSq = null

                        if( snpData.chrNum==dataBestSnp.chrNum ) {
                            rSq = (bestSnpCis==snp) ? 1.0 : LDCalculator.calcLD(snpData,dataBestSnp)
                        }

                        // store LD result
                        ldRsqResults[snp][iso] = new LDResult(snpBestCis:dataBestSnp, snpTrans:snpData, rSq:rSq)

                        // log LD result
                        writer.println("${snp}\t${snpData.locus}\t${bestSnpCis}\t${dataBestSnp.locus}\t${iso}\t${(rSq==null) ? 'NA' : rSq}")
                    } 
                    catch(ex){
                        Logger.getLogger(EqtlStatistics.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                else{
                    System.err.println("Error: best snp eqtl-cis for isoform ${iso} not found.")
                }
            }
        }
        
        writer.close()
    }
    
    /**
     *
     */
    public void generateCorrFileExtended(workDir) {
        def suff = '.trans.filter.ext'
        def fileName = workDir+CORREL_FILE_UNFILT+suff
        
        def writer = new PrintWriter(fileName)
        writer.println( "snp\tisoform" + CORR_FILE_EXTCOLS + CORR_FILE_LDCOLS )
        
        transEqtls.each{ snpId, list->
            SNPData snp = snpData[snpId]
            
            list.each{ isoId->
                def iso = isoforms[isoId]    
                String locus = "${iso.chr}:${iso.start}-${iso.end}"
                BLASTResult blast = ((snp.eqtl) ? blastResults[blastResKey(snp.eqtl,isoId)] : null)
                String blastIdent = ((blast) ? ((blast.alignments) ? blast.averageIdent() : 'NC') : 'Error')
                String blastScore = ((blast) ? ((blast.alignments) ? blast.totalScore() : 'NC') : 'Error')
                LDResult ldRes = ldRsqResults[snpId][isoId] // LD info

                //snp_locus, isoform_locus, iso_eqtl_cis, blast_ident, blast_score, best_snp_iso_cis, best_snp_cis_locus, ld_rsquare
                String extInfoLine = "\t${snp.locus}\t${locus}\t${snp.eqtl ?: 'NA'}\t${blastIdent}\t${blastScore}" +
                    "\t${ldRes?.snpBestCis?.id}\t${ldRes?.snpBestCis?.locus}\t${ldRes?.rSq ?: 'NA'}"
                
                writer.println("${snpId}\t${isoId}" + extInfoLine)
            }
        }
        
        writer.close()
    }
    
}

