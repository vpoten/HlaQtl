/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import weka.core.converters.ArffLoader
import weka.core.Instances
import weka.core.Instance
import weka.core.Attribute
import org.ngsutils.stats.BiNGO.BingoAlgorithm
import org.ngsutils.ontology.GOManager

/**
 *
 * @author victor
 */
class PathwayPatterns {
    
    static RELATION_GENES = 'Genes'
    static RELATION_GO = 'Gene Ontology'
    static OUT_ARFF_GENES = 'patterns_genes.arff'
    static OUT_ARFF_GO = 'patterns_go.arff'
    static final double REPORT_INCR = 10.0
    
    /**
     *
     */
    static perform(args) {
        //get args
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        def pairsFile = Main.readOption(args, Main.OPT_PAIRS)
        def subjsFile = Main.readOption(args, Main.OPT_SUBJECTS)
        def taxonomyId = Main.readOption(args, Main.OPT_TAXID)
        
        if( !pairsFile || !subjsFile || !taxonomyId ){
            System.err.println( "Usage: ${Main.USA_PATH_PATTERNS}")
            return null
        }
        
        
        // load SNP-gene pairs file (see SnpGeneRelation.groovy class)
        def snpGeneMap = [:]
        
        new File(pairsFile).splitEachLine("\\s"){ toks->
            if( snpGeneMap[toks[0]]==null )
                snpGeneMap[toks[0]] = [] as Set
                
            snpGeneMap[toks[0]] << toks[1]
        }
        
        
        // load subjects arff file (see ped_to_arff.groovy script)
        ArffLoader loader = loadDataset(subjsFile)
        Instances data = loader.getDataSet()
        data.setClassIndex(data.numAttributes() - 1)
        
        // check attributes
        for(int i=0; i<data.numAttributes()-1; i++) {
            Attribute att = data.attribute(i)
            assert (att.isNominal()), att.name()+' is not nominal'
            assert (att.numValues()==2), att.name()+' has more than 2 values'
        }
        
        def namespace = GOManager.BP
        GOEnrichment enrichment = new GOEnrichment(taxonomyId:taxonomyId)
        enrichment.loadSemanticData(outdir)
        
        
        def instGOTerms = [] //Selected GO terms per instance
        def instGenes = [] //Selected Genes per instance
        def instClass = []
        double nextReport = REPORT_INCR

        // for each instance perform a GO enrichment analysis
        for(int i=0; i<data.numInstances(); i++) {
            Instance instance = data.instance(i)
            def currentSnps = [] as Set
            
            // filter SNPs
            for(int j=0; j<instance.numAttributes()-1; j++) {
                if( instance.value(j)==1.0 ){
                    currentSnps << data.attribute(j).name()
                }
            }
            
            // get class
            instClass << instance.classValue()
            
            def currentGenes = currentSnps.collect{snpGeneMap[it]}.flatten().findAll{it!=null} as Set
            instGenes <<  currentGenes
            
            // call BiNGO for current instance
            BingoAlgorithm bingo = enrichment.calcBingo(currentGenes)
            
            //filter GO terms by namespace and add subset to instGOterms
            instGOTerms << enrichment.filterByNamespace(bingo, namespace)
            
            // print completion status
            double complete = (i/(double)data.numInstances())*100.0
            if( complete >= nextReport ) {
                println "GO enrichment: ${nextReport}% completed. ${new Date()}"
                nextReport += REPORT_INCR
            }
        }//end of GO enrichment
        
        // write gene patterns
        writeArffFile(outdir+OUT_ARFF_GENES, RELATION_GENES, instGenes, instClass)
        
        // write pathway patterns
        writeArffFile(outdir+OUT_ARFF_GO, RELATION_GO+" ${namespace}", instGOTerms, instClass)
        
    }
    
    /**
     *
     */
    static private def loadDataset(arffFile){
        ArffLoader loader = new ArffLoader()
        loader.setFile(new File(arffFile))
        return loader
    }

    
    /**
     *
     */
    static private def writeArffFile(arffFile, relationName, listFeat, instClass){
        def writer = new BufferedWriter(new FileWriter(arffFile))
        
        writer.writeLine("@RELATION '${relationName}'")
        writer.writeLine('')
        
        def attributes = [] as TreeSet
        listFeat.each{ attributes.addAll(it) }
        
        attributes.each{ writer.writeLine("@ATTRIBUTE '${it}' {0,1}") }
        writer.writeLine("@ATTRIBUTE phenotype {1,2}")
        writer.writeLine('')
        writer.writeLine('@DATA')
        
        listFeat.eachWithIndex{ featSet, i->
            attributes.each{ writer.write( (it in featSet) ? '1,' : '0,' ) }
            writer.writeLine( (instClass[i]==0.0) ? '1' : '2' )
        }
        
        writer.close()
    }
    
}

