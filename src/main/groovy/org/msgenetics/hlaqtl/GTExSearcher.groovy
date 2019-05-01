/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table

import org.ngsutils.Utils
import org.ngsutils.variation.SNPData

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
    
    /** size of region around query SNPs */
    int snpRegionSize = 10000000
    
    /** GTEx eqtl data */
    String gtexDir = null
    
    /** working directory */
    String workDir = null
    
    /** 1000 genomes (vcf + tbi) directory */
    String genomesDir = null
    
    /** GTEx tissues where to filter the eqtls (all by default) */
    List<String> tissues = null
    
    /** eqtl p-value threshold to filter best eqtls */
    double eqtlThr = 0.05d
    
    /** LD threshold to filter results */
    double ldThr = 0.5d
    
    /** keep cache of tped files */
    boolean useCache = true
    
    /** list of subjects to use */
    List<String> subjects
    
    /** Number of threads to use in LD calculation */
    int nThreads = 4
    
    
    /**
     *
     */
    def perform() {
        // Load snp info from ensembl rest api
        EnsemblRestClient client = new EnsemblRestClient(ENSEMBL_REST_API, 15, 200)
        List<SNPData> snpQuery = client.getSnps(queryIds, 'human')
        
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
            Collections.sort(schrRegions[chr], [compare: {a, b -> a.start <=> b.start}] as Comparator)
        }
        
        // Associate eqtls with regions
        chrRegions.each { chr, regions ->
            regions.each { region ->
                Table result = GTExEqtl.filterByRegion(bestEqtls, chr, region.start, region.end)
                regions['eqtls'] = result
            }
        }
        
        // Obtain tped files from 1000genomes vcfs
        chrRegions.each { chr, regions ->
            // TODO use the complete chromosome
            int start = regions.min{ it.start }
            int end = regions.max{ it.end }
            def locusStr = "${chr}:${start}-${end}"
            def vcfFile = SNPManager.S3_VCF_FILE.replace('{chr}', chr)
            def groups = []
            SNPManager.loadSNPData(subjects, workDir, vcfFile, groups, locusStr, true, null)
            // TODO check if the generated .tped file should be renamed
        }
        
        //  Load SNPData from tped files for all snps: query + eqtl
        def snpsData = [:] as TreeMap
        
        chrRegions.each { chr, regions ->
            tpedFile = "${SNPManager._1000G_PED}.tped"
            snps = [] as TreeSet
            
            regions.each { region->
                // add the snp that leads the region and all the snps in associated eqtls
                snps << region.snp.id
                Table eqtls = region.eqtls
                def regionSnps = eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13').asSet().findAll{ it!=null }
                region['region_snps'] = regionSnps
                snps += regionSnps
            }
            
            SNPData.createFromTped(tpedFile, snps, snpsData)
        }
        
        // LD calculation between query snps and eqtl snps in region
        regionLDCalc(chrRegions, snpsData)
        
        // TODO final result
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
                            Double rSq = calcLD(data1, data2)
                            if( rSq!=null ) {
                                region['ld_results'][snpId] = rSq
                            }
                        }
                    }
                }
            }
        }// end lDThreadCalc closure
        
        // run lDThreadCalc in threads
        def closures = (1..nThreads).collect{ (Closure) {lDThreadCalc(it - 1, nThreads)} }
        Utils.runClosures(closures, nThreads )
    }
}

