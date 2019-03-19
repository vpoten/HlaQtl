/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.ngsutils.variation.SNPData
import org.ngsutils.annotation.genome.*
import org.ngsutils.stats.FisherExact

/**
 *
 */
class FisherTestResult {
    String name
    String grp1
    String grp2
    int a
    int b
    int c
    int d
    Double pvalue
    
    String toString(){
        "${name}\t${grp1}\t${grp2}\t${a}\t${b}\t${c}\t${d}\t${pvalue}"
    }
}

/**
 *
 * @author victor
 */
class FeatureSnpsEnrichment {
	
    static final String FISHER_TEST_OUT = 'fisher_tests.out.txt'
    
    /**
     *
     */
    static perform(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String snpsFile = Main.readOption(args, Main.OPT_SNPS)
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        
        println "output: ${outdir}"
        println "snps groups file: ${snpsFile}"
        
        def snpsGroups = [:] as TreeMap //key=group, value=Set of SnpData
        def snpsData = [:] as TreeMap //key=snpId, value=SnpData object
        def taxId = '9606'
        
        // parse snps log extended file
        def reader = new File(snpsFile).newReader()
        reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            def groups = (toks[7].split(',') as Set)
            
            def snp = new SNPData( 
                id:toks[0], chr:toks[1], position:((toks[2] as Integer)-1),
                minor:toks[4], major:toks[5], maf:toks[6] as Double )
            
            //add snp to snpData
            snpsData[toks[0]] = snp
            
            //add snp to each group where it belongs
            groups.each{ grp->
                def set = snpsGroups[grp]
                
                if(set==null){
                    set = [] as TreeSet
                    snpsGroups[grp] = set
                }
                
                set << snp
            }
        }
        
        reader.close()
        
        //create Fisher exact test object
        FisherExact ftest = new FisherExact( snpsGroups.collect{it.value.size()*2}.max() )
        
        //closures section
        def createCounts = { feats->
            def map = [:] as TreeMap
            
            feats.each{ feat->
                def mapCount = [:] as TreeMap
                snpsGroups.keySet().each{mapCount[it] = 0}
                map[feat] = mapCount
            }
            
            return map
        }
        
        ////
        // computes fisher exact test using 2 × 2 contingency table:
        //
        // feature present? | yes | no
        // -----------------------------
        // snp_group_1      | a   | b
        // -----------------------------
        // snp_group_i      | c   | d
        //
        def fisherTests = { mapCounts->
            def results = []
            def groups = snpsGroups.keySet() as List
            
            mapCounts.each{ name, map->
                (1..groups.size()-1).each{
                    int a = map[groups[0]]
                    int b = snpsGroups[groups[0]].size()-a
                    int c = map[groups[it]]
                    int d = snpsGroups[groups[it]].size()-c
                    
                    Double pval = ftest.getCumlativeP(a, b, c, d)
                    
                    results << new FisherTestResult(name:name, grp1:groups[0],
                        grp2:groups[it], a:a, b:b, c:c, d:d, pvalue:pval)
                }
            }
            
            return results
        }
        
        def doCounts = { featClos, countObj->
            snpsData.each{ id, snp->
                def feats = featClos(snp)

                feats.each{ feat->
                    def mapCount = countObj[feat]
                    snpsGroups.each{ grp, set->
                        if( snp in set ){
                            mapCount[grp] = (mapCount[grp]+1)
                        }
                    }
                }
            }
        }
        //end closures
        
        // 1 - Regulation tests
        FeatureIndex regulIndex = 
            AnnotationFactory.createIndex(outdir, taxId, AnnotationFactory.ANN_UCSC_REGUL)
            
        def regulFeatsCounts = createCounts(regulIndex.names.keySet())
        
        doCounts( {snp-> regulIndex.getFeatsByPos("${snp.locus}-${snp.position+1}").collect{it.name}}, 
            regulFeatsCounts )
        
        def resRegul = fisherTests(regulFeatsCounts)
            
        // 2 - Variation tests
        FeatureIndex varIndex = 
            AnnotationFactory.createIndex(outdir, taxId, AnnotationFactory.ANN_UCSC_SNPS138, snpsData.keySet())
            
        def varFeatsCounts = createCounts(AnnotationTrack.SNP_FUNC.keySet())
        
        doCounts( {snp-> varIndex.getFeatsByName(snp.id).collect{it.terms}.flatten()}, 
            varFeatsCounts )
        
        def resVar = fisherTests(varFeatsCounts)
        
        // End - write output table to file
        def outFile = FISHER_TEST_OUT
        def writer = new PrintWriter(outdir+outFile)
        writer.println("name\tgrp1\tgrp2\ta\tb\tc\td\tpvalue\ttype") 
        
        resRegul.each{ writer.println(it.toString()+'\tregulation') }
        resVar.each{ writer.println(it.toString()+'\tvariation') }
        
        writer.close()
        
    }
    
    /**
     *
     */
    static def normalizeChr(taxId, chr){
        return (chr=='23') ? 'X' : chr
    }
    
    /**
     *
     */
    static addFeaturesToSnp(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        String snpsFileStr = Main.readOption(args, Main.OPT_SNPS)
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        
        println "output: ${outdir}"
        println "snps log file: ${snpsFileStr}"
        
        File snpsFile = new File(snpsFileStr)
        def taxId = '9606'
        def snpsData = [:] as TreeMap //key=snpId, value=SnpData object
        
        // parse snps log file
        def reader = snpsFile.newReader()
        def header = reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            def snp = new SNPData( 
                id:toks[0], chr:normalizeChr(taxId, toks[1]), position:((toks[2] as Integer)-1),
                minor:toks[4], major:toks[5], maf:toks[6] as Double )
            
            //add snp to snpData
            snpsData[toks[0]] = snp
        }
        
        reader.close()
        
        // Load regulation and variation
        FeatureIndex regulIndex = 
            AnnotationFactory.createIndex(outdir, taxId, AnnotationFactory.ANN_UCSC_REGUL)
        FeatureIndex varIndex = 
            AnnotationFactory.createIndex(outdir, taxId, AnnotationFactory.ANN_UCSC_SNPS138, snpsData.keySet())
            
        // print output file
        String outFile = snpsFile.name + '.var_reg_added'
        def writer = new PrintWriter(outdir + outFile)
        writer.println(header+"\tregulation\tvariation")
        
        reader = snpsFile.newReader()
        reader.readLine() //skip header
        
        reader.eachLine{ line->
            def toks = line.split("\t",2)
            def snp = snpsData[toks[0]]
            
            def regul = regulIndex.getFeatsByPos("${snp.locus}-${snp.position+1}").collect{it.name} as TreeSet
            def vari = varIndex.getFeatsByName(snp.id).collect{it.terms}.flatten() as TreeSet
            
            writer.println(line+"\t${regul.sum{it+','}}\t${vari.sum{it+','}}")
        }
        
        writer.close()
        reader.close()
    }
    
}

