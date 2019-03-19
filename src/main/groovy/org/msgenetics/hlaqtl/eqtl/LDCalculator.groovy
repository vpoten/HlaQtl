/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.eqtl.LDResult
import org.ngsutils.Utils
import org.ngsutils.variation.HaplotypeFreq
import org.ngsutils.variation.SNPData

/**
 * Linkage desequilibrium standalone calculation
 * 
 * @author victor
 */
class LDCalculator {
    
    static final double ZERO = 1e-12
    
    /**
     * standalone calculation calling plink
     */
    static Map perform(locusStr, snpsQuery, snpsSet, subjects, outdir) {
        def locus = ((locusStr!=null) ? (locusStr =~ Utils.locusRegex) : null)
        def chr = (locusStr!=null) ? locus[0][1] : null
        def start = (locusStr!=null) ? (locus[0][2] as Integer) : null
        def end = (locusStr!=null) ? (locus[0][3] as Integer) : null
        
        // get vcf from 1000genomes ftp using tabix
        def vcfFile = SNPManager._1000G_FTP+SNPManager.S3_VCF_FILE
        def tbiFile = SNPManager._1000G_FTP+SNPManager.S3_VCF_TBI
        vcfFile = vcfFile.replace('{chr}',chr)
        tbiFile = tbiFile.replace('{chr}',chr)
        
        def tmpDir = outdir+'tmp/'
        Utils.createDir(tmpDir)
        String outPref = outdir+SNPManager._1000G_PED+"_${chr}_${start}_${end}"
        def tpedFile = outPref+'.tped'
        
        if( new File(tpedFile).exists() ){
            println "LDCalculator: ${tpedFile} already exists, skip download from 1000genomes."
        }
        else{
            // get file from 1000 genomes using tabix and filter it with vcftools
            SNPManager.vcfToTped(vcfFile, tbiFile, chr, start, end, outPref, 
                false, subjects, tmpDir)
        }
        
        def allSnps = snpsSet
        
        if( !snpsSet ){
            // get snps in tped file
            allSnps = [] as TreeSet
            new File(tpedFile).eachLine{ line->
                def toks = line.split("\t",3)
                allSnps << toks[1]
            }
        }
        
        def snpsData = SNPData.createFromTped(tpedFile, allSnps)
            
        def results = [:] as TreeMap
        def writerLock = new Object()
        
        snpsQuery.each{ id-> //init results map
            results[id] = [:] as TreeMap
        }
        
        //closure for LD calculation (threaded)
        def lDThreadCalc = { thread, totalThrs ->
            def thrLdRsqResults = [:] as TreeMap
            
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
            
            synchronized(writerLock) {
                // add thread results to instance results
                thrLdRsqResults.each{ snp1, map-> map.each{ snp2, res-> results[snp1][snp2] = res } }
            }
        }//end lDThreadCalc closure
        
        //run lDThreadCalc in threads
        int nthreads = 4
        Utils.runClosures(
            [{lDThreadCalc(0, nthreads)}, {lDThreadCalc(1, nthreads)},
            {lDThreadCalc(2, nthreads)}, {lDThreadCalc(3, nthreads)}], 
            nthreads )
        
        return results
    }
	
    /**
     * standalone LD calculation
     */ 
    static public Double calcLD(SNPData data1, SNPData data2) {
        
        if( data1.id==data2.id ) {
            return 1.0d
        }
        
        if( data1.chrNum!=data2.chrNum ){
            //cannot calculate LD
            throw new RuntimeException("Cannot calculate LD. SNPs in different chromosomes")
        }
        
        data1.alleles //build alleles
        data2.alleles
         
        double fA = data1.freqA1
        double fB = data2.freqA1
        double fa = data1.freqA2
        double fb = data2.freqA2
        
        def hapFrq = new HaplotypeFreq(data1, data2)
        hapFrq.eps = 1e-3
        hapFrq.maxIter = 500000
        hapFrq.em()
        
        if( hapFrq.haploFreq==null ){
            System.err.println("LDCalculator: ${data1.id}:${data1.chr} ${data2.id}:${data2.chr} not convergence")
            return -1.0d //error
        }
        
        def freqs = hapFrq.haploFreq
        
        ///double D = freqs[0] - ((freqs[0]+freqs[1])*(freqs[0]+freqs[2]))
        double D = freqs[0] - fA*fB
        
//        if( calculateDp ) {
//          double dmax1 = (D > 0) ? fA * fb : fA * fB;
//          double dmax2 = (D > 0) ? fa * fB : fa * fb;
//          double dmax = (dmax1 < dmax2) ? dmax1 : dmax2;
//          if( dmax<ZERO ){ return -1.0d }
//          return D / dmax;
//        }

        double denom = fA * fa * fB * fb;

        if( denom<ZERO ){ return -1.0d }
        return (D*D) / denom;
    }
    
    
}

