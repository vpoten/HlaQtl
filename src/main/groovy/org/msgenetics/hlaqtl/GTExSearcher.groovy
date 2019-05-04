/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table

import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.IntColumn
import tech.tablesaw.api.DoubleColumn

import org.msgenetics.hlaqtl.eqtl.LDCalculator
import org.msgenetics.hlaqtl.eqtl.SNPManager
import org.msgenetics.hlaqtl.eqtl.GTExEqtl


/**
 *
 * @author victor
 */
class GTExSearcher {
    
    static final String ENSEMBL_REST_API = 'http://grch37.rest.ensembl.org'
    
    /** list of query snps (rs ids) */
    List<String> queryIds = []
    
    /** Size of region around query SNPs */
    int snpRegionSize = 10000000
    
    /** GTEx eqtl data */
    String gtexDir = null
    
    /** Working directory */
    String workDir = null
    
    /** 1000 genomes (vcf + tbi) directory */
    String genomesDir = null
    
    /** GTEx tissues where to filter the eqtls (all by default) */
    List<String> tissues = null
    
    /** eqtl p-value threshold to filter best eqtls */
    double eqtlThr = 0.05d
    
    /** LD threshold to filter results */
    double ldThr = 0.5d
    
    /** Keep cache of tped files */
    boolean useCache = true
    
    /** List of subjects to use */
    List<String> subjects
    
    /** Number of threads to use in LD calculation */
    int nThreads = 4
    
    /** Allowed chromosomes */
    Set<String> chrAllowed = null;
    
    /**
     *
     */
    def perform() {
        // Load snp info from ensembl rest api
        EnsemblRestClient client = new EnsemblRestClient(ENSEMBL_REST_API, 15, 200)
        List<SNPData> snpQuery = client.getSnps(queryIds, 'human')
        
        // apply allowed chromosomes restriction
        if (chrAllowed != null) {
            snpQuery = snpQuery.findAll{it.chr in chrAllowed}
        }
        
        // Load best eqtls from GTEx data
        Table bestEqtls = GTExEqtl.getBestEqtlsAllTissues(gtexDir, eqtlThr);
        
        // Build regions around snps in query list
        def chrRegions = [:]
        
        snpQuery.each { snp ->
                if ( !(snp.chrNum in chrRegions) ) {
                    chrRegions[snp.chrNum] = []
                }
                def start = (snp.position - snpRegionSize) < 1 ? 1 : snp.position - snpRegionSize
                chrRegions[snp.chrNum] << ['start': start, 'end': snp.position + snpRegionSize, 'snp': snp]
            }
        
        // sort chr regions by start position
        for(String chr : chrRegions.keySet()) {
            Collections.sort(chrRegions[chr], [compare: {a, b -> a.start <=> b.start}] as Comparator)
        }
        
        // Associate eqtls with regions
        chrRegions.each { chr, regions ->
            regions.each { region ->
                Table result = GTExEqtl.filterByRegion(bestEqtls, chr, region.start, region.end)
                region['eqtls'] = result
            }
        }
        
        if(useCache == false) {
            // TODO implement check for cache use
            // Obtain tped files from 1000genomes vcfs
            chrRegions.each { chr, regions ->
                // TODO use the complete chromosome
                int start = regions.min{ it.start }.start
                int end = regions.max{ it.end }.end
                def locusStr = "${chr}:${start}-${end}"
                def vcfFile = new File(genomesDir , SNPManager.S3_VCF_FILE.replace('{chr}', "chr${chr}")).absolutePath
                def chrDir = buildChrDir(chr) + '/'
                // generate tped files in a separate directory for each chromosome
                SNPManager.loadSNPData(subjects, chrDir, vcfFile, [], locusStr, true, null)
            }
        }
        
        //  Load SNPData from tped files for all snps: query + eqtl
        def snpsData = [:] as TreeMap
        
        chrRegions.each { chr, regions ->
            def tpedFile = new File(buildChrDir(chr), "${SNPManager._1000G_PED}.tped").absolutePath
            def snps = [] as TreeSet
            
            regions.each { region->
                // add the snp that leads the region and all the snps in associated eqtls
                snps << region.snp.id
                Table eqtls = region.eqtls
                StringColumn  regionSnps = (StringColumn) eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13').asSet().findAll{ it!=null }
                region['region_snps'] = regionSnps
                snps += regionSnps
            }
            
            SNPData.createFromTped(tpedFile, snps, snpsData)
        }
        
        // LD calculation between query snps and eqtl snps in region
        regionLDCalc(chrRegions, snpsData)
        
        // Final result
        Table finalTable = Table.create("GTEx_eqtl_filter_LD")

        chrRegions.each { chr, regions ->
            regions.each { region->
                // get snpIds with rSq > threshold
                def ldResultsPass = region['ld_results'].findAll({ snpId, rSq -> rSq > ldThr})
                
                // create the final result table for the region: eqtl table + extra columns:
                // extra columns: region_snp, region_snp_pos, ld_rsq
                Table eqtls = region.eqtls
                StringColumn regionSnps = (StringColumn) eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13')
                Table result = eqtls.filter(regionSnps.isIn(ldResultsPass.keySet()))
                
                StringColumn resultSnps = (StringColumn) result.stringColumn('rs_id_dbSNP147_GRCh37p13')
                // create extra columns
                def snpCol = StringColumn.create('region_snp', (0..result.rowCount()-1).collect{region.snp.id})
                def snpPosCol = IntColumn.create('region_snp_pos', (0..result.rowCount()-1).collect{region.snp.position})
                def ldRsqCol = DoubleColumn.create('ld_rsq', resultSnps.asList().collect{region['ld_results'][it]})
                
                // add the columns to the result table
                result.addColumns(snpCol, snpPosCol, ldRsqCol)
                
                // append to final table
                finalTable.append(result)
            }
        }
        
        finalTable.write().csv(new File(workDir, "${finalTable.name()}_${new Date().format('yyyyMMddHHmmss')}.csv"))
    }
    
    /**
     *
     */
    private def regionLDCalc(chrRegions, snpsData) {
        // collect all regions
        def allRegions = chrRegions.collect({chr, regions -> regions}).flatten()
        
        // closure for LD calculation (threaded)
        def lDThreadCalc = { thread, totalThrs ->
            // assign index to region to divide jobs among threads
            allRegions.eachWithIndex{ region, i-> 
                if ((i % totalThrs) == thread) {
                    region['ld_results'] = [:] as TreeMap
                    def data1 = region.snp
                    def current = region.snp.id

                    region['region_snps'].each { snpId ->
                        if( current != snpId ) {
                            def data2 = snpsData[snpId]
                            Double rSq = LDCalculator.calcLD(data1, data2)
                            if( rSq!=null ) {
                                region['ld_results'][snpId] = rSq
                            }
                        }
                    }
                }
            }
        }// end lDThreadCalc closure
        
        // run lDThreadCalc in threads
        nThreads = 4
        
        def closures = [
            {lDThreadCalc(0, nThreads)},
            {lDThreadCalc(1, nThreads)},
            {lDThreadCalc(2, nThreads)},
            {lDThreadCalc(3, nThreads)}
        ]
        
        Utils.runClosures(closures, nThreads )
    }
    
    /**
     * Set subjects to use in calculations
     */
    def setSubjects(populations) {
        subjects = IGSRSamples.create().getSubjects(populations)
    }
    
    /**
     * Load query snp ids from text file
     */
    def loadQuerySnpsFromFile(file) {
        def reader = new File(file).newReader()
        queryIds.clear()
        
        reader.eachLine { line ->
            queryIds << line.trim()
        }
        
        reader.close()
    }
    
    /**
     * Build a directory for a chromosome within workDir and returns the created path
     */
    def buildChrDir(String chr) {
        if( !chr.startsWith('chr') ){
            chr = 'chr' + chr
        }
        def chrDir = new File(workDir, chr)
        String chrPath = chrDir.absolutePath
        
        if (!chrDir.exists()) {
            Utils.createDir(chrPath)
        }
        
        return chrPath
    }
    
    /**
     * Set allowed chromosomes to use in calculations
     */
    def setChrAllowed(chrs) {
        chrAllowed = chrs as Set
    }
}

