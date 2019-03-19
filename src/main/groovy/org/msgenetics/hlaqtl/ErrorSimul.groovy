/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl


import org.msgenetics.hlaqtl.eqtl.EqtlSearcher as ESe
import org.ngsutils.variation.SNPData
import org.msgenetics.hlaqtl.eqtl.IsoformParser
import org.ngsutils.eqtl.CorrelationCalc
import org.ngsutils.Utils
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

/**
 * Simulation of eQTLs correlation calculation adding typing error to genotypes.
 * Measures the difference between the original correlation value and altered 
 * genotype correlation value, dividing the eQTLs in intervals of correlation 
 * value and MAF.
 * 
 * @author victor
 */
class ErrorSimul {
    
    static double corrStep = 0.1
    static double mafStep = 0.1
    static double typingError = 0.0025
    static double missingRate = 0.0
    static int limitCorr = 10000
    static int nthreads = 2
    
    static ESe searcher = null //EqtlSearcher
    static def subjects //subject list (IDs)
    static def isoExprData = [:] as TreeMap //key=subjectId, value={map with key=iso,value=expresion}
    
    // Constants
    static final String SNP_LOG = 'snps_info.log'
	
    /**
     *
     */
    static perform(args, int p_nthreads) {
        nthreads = p_nthreads
        
        //read input and output dir
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String eqtlCisDir = args.find{ it.startsWith(Main.OPT_EQTL_DIR) }
        String isoFile = Main.readOption(args, Main.OPT_ISOFORMS)
        String limitStr = Main.readOption(args, Main.OPT_VALUE)
        
        def dirs = [outdir, eqtlCisDir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        eqtlCisDir = dirs[1]
        
        println "output: ${outdir}"
        println "input eqtls-cis dir: ${eqtlCisDir}"
        println "isoforms filtered file: ${isoFile}"
        println "typing error rate: ${typingError}"
        
        if( limitStr ){
            limitCorr = limitStr as Integer
            println "Number of eQTLS per interval: ${limitCorr}"
        }
        
        // parse isoforms filtered file
        def isoforms = IsoformParser.parseIsoFilterFile(isoFile)
        
        // search eQTLs
        searcher = new ESe()//new EqtlSearcher
        def isoCorrMap = searcher.searchIsoformsExtended(eqtlCisDir, isoforms)
        
        println "\n${isoforms.size()} isoforms read before eQTLs filter"
        println "${searcher.snpData.size()} SNPs read after eQTLs filter"
        println "${isoCorrMap.size()} isoforms read after eQTLs filter\n"
        
        File snpLog = new File(outdir+SNP_LOG)
        
        if( !snpLog.exists() ) {
            // calc SNPs frequencies
            println "Calculating SNPs frequencies. ${new Date()}\n"
            searcher.calSnpFrequencies()
            genSnpsInfoLog(snpLog, searcher.snpData) 
        }
        else {
            println "Parsing previously generated SNPs frequencies log file.\n"
            parseSnpsInfoLog(snpLog, searcher.snpData)
        }
        
        subjects = SNPData.subjectList( searcher.snpData )
        println "${subjects.size()} subjects read."
        
        println "Parsing isoform expression data.\n"
        parseIsoformTracking(eqtlCisDir, isoforms)
        
        // 1 - assign eQTLs to each correlation interval
        def corrIntervals = (1..(1.0/corrStep)).collect{ [] }

        isoCorrMap.each{ iso, list->
            list.each{ eqtlData->
                int idx = eqtlData[ESe.F_CORR]/corrStep
                corrIntervals[idx] << eqtlData
            }
        }
        
        printIntervals('Correlation', corrIntervals, corrStep)
        
        // 2 - assign eQTLs to each MAF interval
        def mafIntervals = (1..5).collect{ [] }
        
        isoCorrMap.each{ iso, list->
            list.each{ eqtlData->
                def snpdata = searcher.getSnpData(eqtlData[ESe.F_SNP])
                int idx = snpdata.maf/mafStep
                
                if(idx==5){
                    ///System.err.println("SNP freq error: ${snpdata.id},${snpdata.locus},${snpdata.alleles},${snpdata.maf}")
                    idx=4
                }
                
                mafIntervals[idx] << eqtlData
            }
        }
        
        printIntervals('MAF', mafIntervals, mafStep)
        
        // 3 - count missing alleles and calc. missing rates
        int missSnps = SNPData.countMissing( searcher.snpData )
        long totalBases = 2L*subjects.size()*searcher.snpData.size()
        missingRate = missSnps/(double)totalBases
        
        println "Total amount of bases = ${totalBases}"
        println "Missing SNPs = ${missSnps}"
        println "Missing rate = ${missingRate}"

        
        def correlation = { title, intervals->
            intervals.eachWithIndex{ list, i->
                def selectedIdxs = randomIndexSet(limitCorr, list.size())
                calcCorrelation(list, selectedIdxs, outdir, "${title}_${i}")
            }
        }//end closure
        
        ['Corr':corrIntervals, 'MAF':mafIntervals].each{ t, i-> correlation(t,i) }
    }
    
    
    /**
     * print out count of eQTLs per interval
     */
    private static printIntervals(title, intervals, step){
        intervals.eachWithIndex{ list, i->
            println "${title}: interval ${i} (<${(i+1)*step}) size = ${list.size()}"
        }
        
        println "${title}: ${intervals.sum{it.size()}} eQTLs read after filtering"
    }
    
    
    /**
     *
     */
    private static genSnpsInfoLog(File file, snpDataMap){
        def writer = new PrintWriter(file)
        writer.println("id\tchr\tposition\talleles\tminor\tmajor\tmaf")
        
        snpDataMap.each{ id, data->
            writer.println("${id}\t${data.chr}\t${data.position}\t${data.alleles}\t${data.minor}\t${data.major}\t${data.maf}")
        }
        
        writer.close()
    }
    
    
    /**
     *
     */
    private static parseSnpsInfoLog(File file, snpDataMap){
        def reader = file.newReader()
        reader.readLine()//skip header
        
        reader.splitEachLine("\t"){ toks->
            def data = snpDataMap[toks[0]]
            data.alleles = toks[3]
            data.minor = toks[4]
            data.major = toks[5]
            data.maf = toks[6] as Double
        }
        
        reader.close()
    }
    
    /**
     *
     */
    private static genCorrOutFile(outdir, title){
        "${outdir}${Main.CORREL_FILE}.${title}"
    }
    
    /**
     *
     */
    private static calcCorrelation(listEqtls, selectedIdxs, outdir, title){
        def snpSet = [] as TreeSet
        def isoSet  = [] as TreeSet
        def corrSet = [] as TreeSet
        
        selectedIdxs.each{
            def eqtl = listEqtls[it]
            snpSet << eqtl[ESe.F_SNP]
            isoSet << eqtl[ESe.F_ISO]
            corrSet << CorrelationCalc.buildKey(eqtl[ESe.F_SNP], eqtl[ESe.F_ISO])
        }
        
        def corrCalc = new CorrelationCalc( subjects: subjects, 
            snps:snpSet, isoforms:isoSet, threads:nthreads )
        
        //populate CorrelationCalc genotypes and expression data
        subjects.eachWithIndex{ subjId, idx->
            def mapGeno = [:] as TreeMap
            def mapExpr = [:] as TreeMap
            
            //get genotypes
            snpSet.each{ snpId->
                SNPData snpData = searcher.getSnpData(snpId)
                def alleles = snpData.getSubjectAlleles(idx)
                Double value = snpData.mutateEncode(alleles, typingError, missingRate)
                if( value!=null )
                    mapGeno[snpId] = value
            }
            
            //get expression
            def subjIsoforms = isoExprData[subjId]
            isoSet.each{ isoId->
                Double expr = subjIsoforms[isoId]
                mapExpr[isoId] = ((expr!=null) ? expr : 0.0) as Double
            }
            
            corrCalc.genotypes[subjId] = mapGeno
            corrCalc.expression[subjId] = mapExpr
        }
       
        println "\n${title}: Perform correlation in ${corrCalc.threads} threads. ${new Date()}"
        println "${title}: ${corrCalc.subjects.size()} subjects"
        println "${title}: ${corrCalc.snps.size()} SNPs"
        println "${title}: ${corrCalc.isoforms.size()} isoforms"
        println "${title}: ${corrSet.size()} eQTLs"
        corrCalc.calcSparse(corrSet)
        
        //generate text results
        String outFile = genCorrOutFile(outdir,title)
        corrCalc.printResults( new PrintWriter(outFile) )
        
        //generate compare file and calculate summary statistics
        String compFile = outFile+'.comp'
        def writer = new PrintWriter(compFile)
        def summary = new SummaryStatistics()
        
        def reader = new File(outFile).newReader()
        reader.readLine()
        writer.println("snp\tiso\tcorr1\tcorr2\tdiff")
        
        reader.splitEachLine("\t"){ toks->
            def snp = toks[ESe.F_SNP]
            def iso = toks[ESe.F_ISO]
            def corr2 = toks[ESe.F_CORR] as Double
            def idx = selectedIdxs.find{ (listEqtls[it][ESe.F_SNP]==snp && listEqtls[it][ESe.F_ISO]==iso) }
            double corr1 = listEqtls[idx][ESe.F_CORR]
            summary.addValue(corr1-corr2)
            
            writer.println("${snp}\t${iso}\t${corr1}\t${corr2}\t${corr1-corr2}")
        }
        
        reader.close()
        writer.close()
        
        println "Summary statistics for ${title}:"
        println summary.toString()
    }
    
    
    /**
     * parses isoform tracking files and fills isoExprData map
     * 
     */
    static parseIsoformTracking(dir, selectedIso){
        subjects.each{ isoExprData[it] = [:] as TreeMap }
        
        new File(dir).eachDir{ fDir->
            File mergTrackFile = Utils.getFileIfCompress(fDir.absolutePath+'/'+Main.MERG_TRACK_FILE)

            if( mergTrackFile ){
                def reader = Utils.createReader(mergTrackFile)
        
                reader.splitEachLine("\\s"){ toks->
                    //file fields: subj iso iso_merg locus fpkm
                    if( (toks[2] in selectedIso) && isoExprData.containsKey(toks[0]) ){
                        isoExprData[toks[0]][toks[2]] = toks[4] as Double
                    }
                }

                reader.close()
            }
        }
    }
    
    /**
     *  @param size : size of returned indexes set (< limit)
     *  @param limit : max index to include in set
     *  @return a Set of randomly generated integers between 0 and limit-1
     */
    static Set randomIndexSet( int size, int limit) {
        def set = [] as TreeSet
        
        if(limit==0)
            return set
            
        if( size>=limit ){
            (0..limit-1).each{ set << it }
            return set
        }
            
        while( set.size()<size ){
            set << (int)Math.floor(Math.random()*limit)
        }
        return set
    }
    
}
