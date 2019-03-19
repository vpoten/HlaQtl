/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.msgenetics.hlaqtl.eqtl.EqtlSearcher
import org.ngsutils.semantic.LinkedLifeDataFactory
import org.ngsutils.semantic.query.GeneQueryUtils

/**
 * Class that calculates the relation between a list of SNPs and expressed genes, 
 * using previously calculated eQTLs. As a result of the process a list of pairs 
 * SNP-gene is obtained. The result of this process is used as input data for
 * PathwayPatterns algorithm.
 * 
 * @author victor
 */
class SnpGeneRelation {
	
    
    public static final String SNP_GENE_REL_FILE = 'snp_gene_rel.txt'
    public static final MISSING_CODES = ['NA','?','.','-','0']
    
    
    /**
     *
     */
    static perform(args) {
        //read input and output dir
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(Main.OPT_INPUT_DIR) }
        
        def dirs = [outdir,inputdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        inputdir = dirs[1]
        
        def snpsFile = Main.readOption(args, Main.OPT_SNPS)
        def pairsFile = Main.readOption(args, Main.OPT_PAIRS)
        def taxonomyId = Main.readOption(args, Main.OPT_TAXID)
        
        if( !snpsFile ){
            System.err.println( "Usage: ${Main.USA_SNP_GENE}")
            return null
        }
        
        println "output: ${outdir}"
        println "input eqtls dir: ${inputdir}"// eqtls results dir
        println "snps file: ${snpsFile}"
        println "taxonomy id: ${taxonomyId}"
        
        if( pairsFile ){
            // known association between snps and genes (external source)
            println "pairs file: ${pairsFile}"
        }
        
        Map trueEqtls = EqtlSearcher.calcTrueEqtls(snpsFile, inputdir, outdir)
        
        def nonEqtls = [] as TreeSet
        trueEqtls.each{ if(!it.value){ nonEqtls<<it.key } }
        
        if( pairsFile && nonEqtls ) {
            println "Assigning genes to ${nonEqtls.size()} snps without any reported eQTL"
            
            // assign pairs to non-eqtl snps using the external source
            new File(pairsFile).splitEachLine("\\s"){ toks->
                if( (toks[0] in nonEqtls) && !(toks[1] in MISSING_CODES) ){
                    trueEqtls[toks[0]] << toks[1]
                }
            }
        }
        
        // load entrez gene semantic data
        println 'Loading semantic data ...'
        def graph = LinkedLifeDataFactory.loadRepository(LinkedLifeDataFactory.LIST_EGENE, [taxonomyId], outdir)
        GeneQueryUtils geneQuery = new GeneQueryUtils(graph)
        
        // generate snp-gene relation output file
        println "\nEnding: generating ${outdir+SNP_GENE_REL_FILE}"
        def writer = new PrintWriter(outdir+SNP_GENE_REL_FILE)
        
        trueEqtls.each{ snp, list->
            def symbols = [] as Set
            
            list.each{ iso->
                def uri = geneQuery.getGeneByName(iso, taxonomyId)
                
                if( uri!=null ){
                    def gene = geneQuery.getGeneSymbol(uri)

                    if( gene!=null && !(gene in symbols) ){
                        writer.println("${snp}\t${gene}") 
                        symbols << gene
                    }
                    else if(gene==null){
                        println "Cannot find gene symbol for ${uri.stringValue()}"
                    }
                }
                else {
                    println "Cannot find gene ID for name ${iso}"
                }
            }
        }
        
        writer.close()
    }
    
}

