/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import org.ngsutils.eqtl.CorrelationCalc
import org.ngsutils.variation.PopAnalysis

/**
 *
 * @author victor
 */
class EqtlSimpleCalc {

    def snps
    def subjGeno = [:] as TreeMap //map with key=subjId and value=array of genotypes (encoded)
    def features // expresed genes or transcripts
    def subjExpr = [:] as TreeMap //map with key=subjId and value=array of expression values
    int nthreads = 2
    
    /**
     * 
     */ 
    static performCalc(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def locus = Main.readOption(args, Main.OPT_LOCUS)
        def groups = Main.getOptGroups(args)
        def exprFile = Main.readOption(args, Main.OPT_EXPRDATA)
        String value = Main.readOption(args, Main.OPT_VALUE)
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        
        outdir = dirs[0]
        
        if( !locus || !exprFile ) {
            System.err.println("${Main.USA_EQTL_SIMPLE}")
            System.exit(1)
        }
       
        println "locus: ${locus}"
        println "groups: ${groups}"
        println "output: ${outdir}"
        println "expression file: ${exprFile}"
        println "mode: ${value}"
        
        boolean _1000gOnly = value.toUpperCase()=='HAPMAP' ? false : true;
        
        if( !_1000gOnly && !groups ){ groups = ['CEU'] }
        
        EqtlSimpleCalc instance = new EqtlSimpleCalc()
        instance.loadExpression(exprFile) //first: load expression from file
        
        
        instance.loadGenotypes(null, outdir, locus, _1000gOnly, groups ?: []);
        
        instance.calcEqtls(outdir)
    }
    
    /**
     * 
     */ 
    def calcEqtls(tpedFile, mergedTrkFile, outdir, boolean adjustGeno) {
        loadMergedTracking(mergedTrkFile)
        loadGenotypesFromTped(tpedFile)
        
        if( adjustGeno ){
            adjustGenotypes()
        }
        
        calcEqtls(outdir)
    }

    /**
     *
     */
    def calcEqtls(outdir) {
        def corrCalc = new CorrelationCalc( subjects:subjExpr.keySet(), 
            snps:snps, isoforms:features, threads:nthreads )
        
        //populate CorrelationCalc genotypes and expression data
        subjExpr.keySet().each{ subjId->
            def mapGeno = [:] as TreeMap
            def mapExpr = [:] as TreeMap
            
            snps.eachWithIndex{ id, i->  mapGeno[id] = subjGeno[subjId][i] }
            features.eachWithIndex{ id, i->  mapExpr[id] = subjExpr[subjId][i] }
            
            corrCalc.genotypes[subjId] = mapGeno
            corrCalc.expression[subjId] = mapExpr
        }
        
        Main.debugCorrData(corrCalc, outdir)
        
        println "${Main.LINE_PRE}Perform correlation in ${corrCalc.threads} threads. ${new Date()}"
        println "${corrCalc.subjects.size()} subjects"
        println "${corrCalc.snps.size()} SNPs"
        println "${corrCalc.isoforms.size()} isoforms"
        corrCalc.calcCorrelations()
        
        //generate text results
        println "${Main.LINE_PRE}Generating correlation results in ${Main.CORREL_FILE}. ${new Date()}"
        corrCalc.printResults( new PrintWriter(outdir+Main.CORREL_FILE) )
        
        // filter correlation results file
        double pvalThr = 0.05
        Main.filterCorrResults(outdir+Main.CORREL_FILE, Main.MIN_OBS, pvalThr)
    }
    
    /**
     *
     */
    def loadGenotypes(vcfFile, outdir, locus, boolean _1000gOnly, groups = []) {
        
        SNPManager snpMng = SNPManager.loadSNPData(subjExpr.keySet(), outdir,
            vcfFile, groups, locus, _1000gOnly, null, true, null)
        
        snps = snpMng.imputedSnps.keySet().collect{it}
        
        subjExpr.keySet().each{ subjId->
            def genoObj = snpMng.genotypes[subjId]
            
            subjGeno[subjId] = snps.collect{ snpId->
                def alleles = genoObj?.snps?.get(snpId)
                (alleles) ? snpMng.imputedSnps[snpId].encode(alleles) : null
            }
        }
    }
    
    /**
     * 
     */
    def loadGenotypesFromTped(tpedFile) {
        List tfamSubjs = SNPData.subjectList(new File(tpedFile))
        
        def subjIdxs = [:] as TreeMap
        tfamSubjs.eachWithIndex{val, i-> 
            def id = val.substring( val.indexOf(':')+1 )
            subjIdxs[id] = i
        }
        
        // get snps ids from tped file
        def reader = Utils.createReader(new File(tpedFile))
        def snpsSet = [] as TreeSet
        
        reader.eachLine{ line->
            def toks = line.split("\\s",5)
            snpsSet << toks[1] //rsId
        }
        
        reader.close()
        
        def snpMap = SNPData.createFromTped(tpedFile, snpsSet)
        snps = snpsSet
        
        subjExpr.keySet().each{ subjId->
            subjGeno[subjId] = snps.collect{ snpId->
                def snpData = snpMap[snpId]
                def idx = subjIdxs[subjId]
                
                (idx!=null) ? snpData.encode( snpData.getSubjectAlleles(idx) ) : null
            }
        }
    }
    
    /**
     * adjust genotypes using PCA (Eigenstrat)
     */
    def adjustGenotypes() {
        PopAnalysis instance = new PopAnalysis(subjExpr.keySet() as List, snps.size(), nthreads)
        instance.numAxes = 2
        
        subjGeno.eachWithIndex{ subj, geno, i->
            geno.eachWithIndex{ val, j->
                instance.setGenotype(i,j,val)
            }
        }
        
        instance.eigenAdjustment()
      
        subjGeno.eachWithIndex{ subj, geno, i->
            (0..geno.size()-1).each{ j->
                def val = instance.getGenotype(i,j)
                geno[j] = (val==PopAnalysis.MISSING) ? null : val
            }
        }
    }
    
    /**
     *
     */
    def loadHapMapGenotypes(outdir, locus, groups) {
        
        SNPManager snpMng = SNPManager.loadHapMapSNPData(subjExpr.keySet(), 
            outdir, groups, locus)
        
        snps = snpMng.imputedSnps.keySet().collect{it}
        
        subjExpr.keySet().each{ subjId->
            def genoObj = snpMng.genotypes[subjId]
            
            subjGeno[subjId] = snps.collect{ snpId->
                def alleles = genoObj?.snps?.get(snpId)
                (alleles) ? snpMng.imputedSnps[snpId].encode(alleles) : null
            }
        }
    }
    
    
    /**
     * loads expression from a tab separated file with the following line format:
     * subj_id expr_1 expr_2 ... expr_n
     */
    def loadExpression(file) {
        def reader = Utils.createReader(new File(file))
        def header = reader.readLine()
        
        def tokArr = header.split("\t")
        // gather genes/isoforms
        features = (1..tokArr.length-1).collect{tokArr[it]}
        
        reader.splitEachLine("\t"){ toks->
            subjExpr[toks[0]] = (1..features.size()).collect{ 
                    def str = toks[it].trim()
                    (!str || str=='NA') ? null : (str as Double)
                }
        }
        
        reader.close()
    }
    
    /**
     * loads expression from merged tracking file (generated by previous eqtl calculation)
     */
    def loadMergedTracking(file) {
        def subjects = [:] as TreeMap
        
        // load subjects map
        def reader = Utils.createReader(file)
        
        reader.splitEachLine("\\s"){ toks->
            //subj iso iso_merg locus fpkm
            subjects[toks[0]] = toks[0]
        }
        
        reader.close()
        
        def listIsoData = IsoformParser.parseSubjectIsoforms(file, subjects, null)
        features = listIsoData[0].isoforms.keySet() as List
        
        
        listIsoData.each{ isoParser->
            subjExpr[isoParser.subject] = features.collect{ isoParser.isoforms[it].fpkm }
        }
    }
    
}

