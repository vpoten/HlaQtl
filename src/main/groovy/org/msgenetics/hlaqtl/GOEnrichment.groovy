/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.ngsutils.Utils
import org.ngsutils.stats.BiNGO.BingoAlgorithm
import org.ngsutils.stats.StatTestParams
import org.ngsutils.semantic.NGSDataResource
import org.ngsutils.semantic.LinkedLifeDataFactory
import org.ngsutils.ontology.GOManager
import org.ngsutils.semantic.query.GeneQueryUtils

/**
 *
 * @author victor
 */
class GOEnrichment {
    
    static final String GO_ENRICH_OUT = 'go_enrichment.out.txt'
    
    GeneQueryUtils geneQuery
    NGSDataResource annotation
    GOManager goManager
    def taxonomyId

    /**
     *
     */
    def loadSemanticData(workdir) {
        def graph = LinkedLifeDataFactory.loadRepository(LinkedLifeDataFactory.LIST_BASIC_GO, [taxonomyId], workdir)
        geneQuery = new GeneQueryUtils(graph)
        
        //prepare annotation
        goManager = new GOManager(graph)
        annotation = NGSDataResource.create(graph)
        annotation.load(taxonomyId)
    }
    
    /**
     *
     */
    BingoAlgorithm calcBingo(Set featureSet) {
        
        //prepare stats params
        StatTestParams statParams = new StatTestParams()
        statParams.annotation = annotation
        statParams.taxonomyId = taxonomyId
        ///statParams.namespace = GOManager.BP
        
        //translate features to genes URI
        def genesAsURIs = 
                featureSet.collect{ geneQuery.getGeneByName(it, taxonomyId)?.stringValue() }.findAll{ it!=null }
            
        statParams.annotation.selectedGenes = genesAsURIs

        // call BiNGO for current set
        BingoAlgorithm enrichment = BingoAlgorithm.performCalculations(statParams)
    
        return enrichment
    }
    
    /**
     *
     */
    def filterByNamespace(BingoAlgorithm bingo, namespace) {
        return bingo.correctionMap.findAll{ goManager.getNamespace(it.key)==namespace }.keySet()
    }
    
    /**
     * Load sets from a pairs file
     */
    static def loadSets(pairsFile, boolean header = false) {
        Map sets = [:]
        
        def reader = Utils.createReader(new File(pairsFile))
        if( header ){ reader.readLine() }
        
        reader.splitEachLine("\t"){ toks->
            //<feature> <group_label>
            if( !sets.containsKey(toks[1]) ){
                sets[toks[1]] = [] as Set
            }
            
            sets[toks[1]] << toks[0]
        }
        
        reader.close()
        
        return sets
    }
 
    /**
     *
     */
    static perform(args) {
        //get args
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        def taxonomyId = Main.readOption(args, Main.OPT_TAXID)
        def pairsFile = Main.readOption(args, Main.OPT_PAIRS)
        
        if( !pairsFile || !taxonomyId ){
            System.err.println( "Usage: ${Main.USA_GO_ENRICH}")
            return null
        }
        
        GOEnrichment enrichment = new GOEnrichment(taxonomyId:taxonomyId)
        enrichment.loadSemanticData(outdir)
        
        // collect sets
        Map sets = loadSets(pairsFile)
        
        //perform enrichment analysis for each set
        Map results = [:]
        
        sets.each{ label, features->
            BingoAlgorithm bingo = enrichment.calcBingo(features)
            results[label] = bingo
        }
        
        // write output file
        def outFile = GO_ENRICH_OUT
        def writer = new PrintWriter(outdir+outFile)
        writer.println("GO_term\tp_val\tnamespace\tset") 
        
        results.each{ label, bingo->
            bingo.correctionMap.each{ term, pval->
                def namespace = enrichment.goManager.getNamespace(term)
                writer.println("${term}\t${pval}\t${namespace}\t${label}") 
            }
        }
        
        writer.close()
        
    }
    
    
}

