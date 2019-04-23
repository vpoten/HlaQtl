/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table

import org.ngsutils.variation.SNPData
import org.msgenetics.hlaqtl.eqtl.GTExEqtl


/**
 *
 * @author victor
 */
class GTExSearcher {
    
    static final String ENSEMBL_REST_API = 'http://grch37.rest.ensembl.org'
    
    // list of query snps (rs ids)
    List<String> queryIds = []
    
    // size of region around query SNPs
    int snpRegionSize = 10000000
    
    // GTEx eqtl data
    String gtexDir = null
    
    // working directory
    String  workDir = null
    
    // GTEx tissues where to filter the eqtls (all by default)
    List<String> tissues = null
    
    // Eqtl p-value threshold to filter best eqtls
    double eqtlThr = 0.05d
    
    // LD threshold to filter resuklts
    double ldThr = 0.5d
    
    
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
                chrRegions[snp.chrNum] << new Tuple(start, snp.position + snpRegionSize, {})
            }
        
        // sort chr regions by start position
        for(String chr : chrRegions.keySet()) {
            Collections.sort(schrRegions[chr], [compare: {a, b -> a[0] <=> b[0]}] as Comparator)
        }
        
        // Associate eqtls with regions
        chrRegions.each { chr, regions ->
            regions.each { region ->
                Table result = GTExEqtl.filterByRegion(bestEqtls, chr, region[0], region[1])
                regions[3]['eqtls'] = result
            }
        }
        
        // TODO Obtain tped files from 1000genomes vcfs
        
        // TODO Load SNPData from tped files for all snps: query + eqtl
        
        // TODO LD calculation between query snps and eqtl snps in region
        
        // TODO final result
    }
}

