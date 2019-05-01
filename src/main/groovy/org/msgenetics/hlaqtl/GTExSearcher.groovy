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
                chrRegions[snp.chrNum] << new Tuple(start, snp.position + snpRegionSize, ['snp': snp])
            }
        
        // sort chr regions by start position
        for(String chr : chrRegions.keySet()) {
            Collections.sort(schrRegions[chr], [compare: {a, b -> a[0] <=> b[0]}] as Comparator)
        }
        
        // Associate eqtls with regions
        chrRegions.each { chr, regions ->
            regions.each { region ->
                Table result = GTExEqtl.filterByRegion(bestEqtls, chr, region[0], region[1])
                regions[2]['eqtls'] = result
            }
        }
        
        // Obtain tped files from 1000genomes vcfs
        chrRegions.each { chr, regions ->
            int start = regions.min{ it[0] }
            int end = regions.max{ it[1] }
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
                snps << region[2].snp.id
                Table eqtls = region[2].eqtls
                def regionSnps = eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13').asSet().findAll{ it!=null }
                region[2]['region_snps'] = regionSnps
                snps += regionSnps
            }
            
            SNPData.createFromTped(tpedFile, snps, snpsData)
        }
        
        // LD calculation between query snps and eqtl snps in region
        results = regionLDCalc(chrRegions, snpsData)
        
        // TODO final result
    }
    
    /**
     *
     */
    private def regionLDCalc(chrRegions, snpsData) {
        def results = [:] as TreeMap
        def writerLock = new Object()
        
        // TODO save results in region
        snpsQuery.each{ id-> //init results map
            results[id] = [:] as TreeMap
        }
        
        // TODO collect all regions
        
        //closure for LD calculation (threaded)
        def lDThreadCalc = { thread, totalThrs ->
            def thrLdRsqResults = [:] as TreeMap
            
            // TODO assign index to region to divide jobs among threads
            snpsQuery.each{ current-> 
                thrLdRsqResults[current] = [:] as TreeMap
                def data1 = snpsData[current]

                allSnps.eachWithIndex{ snpId, i->
                    if( (i%totalThrs)==thread && current!=snpId ){
                        def data2 = snpsData[snpId]
                        Double rSq = calcLD(data1, data2)

                        if( rSq!=null ){ thrLdRsqResults[current][snpId] = rSq }
                    }
                }
            }
            
            // TODO will be not necessary
            synchronized(writerLock) {
                // add thread results to instance results
                thrLdRsqResults.each{ snp1, map-> map.each{ snp2, res-> results[snp1][snp2] = res } }
            }
        }//end lDThreadCalc closure
        
        //run lDThreadCalc in threads
        def closures = (1..nThreads).collect{ (Closure) {lDThreadCalc(it - 1, nThreads)} }
        Utils.runClosures(closures, nThreads )
        
        return results
    }
}

