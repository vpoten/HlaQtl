/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.msgenetics.hlaqtl.eqtl.EqtlStatistics
import org.ngsutils.Utils
import org.ngsutils.EnsemblUtils
import org.ngsutils.semantic.query.GeneQueryUtils
import org.ngsutils.semantic.LinkedLifeDataFactory

/**
 *
 * @author victor
 */
class DownSequences {
    
    
    static perform(args) {
        //read input and output dir
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String isoFile = Main.readOption(args, Main.OPT_VALUE)
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        
        println "output: ${outdir}"
        println "input isoform file: ${isoFile}"
        
        //load semantic data
        println 'Loading semantic entrezgene data ...'
        def taxId = '9606'
        def graph = LinkedLifeDataFactory.loadRepository(LinkedLifeDataFactory.LIST_EGENE, [taxId], outdir)
        def geneQuery = new GeneQueryUtils(graph)
        
        
        //process isoform file
        def reader = new File(isoFile).newReader()
        int isoIdx = 0
        
        File seqDir = new File(outdir+EqtlStatistics.SEQS_DIR)
        Utils.createDir(seqDir.absolutePath)
        
        reader.splitEachLine("\\s"){ toks->
            def isoId = toks[isoIdx]
            def geneUri = geneQuery.getGeneByName(isoId, taxId)
                
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
        
        reader.close()
        
    }
    
    
}

