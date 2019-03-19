/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import org.ngsutils.eqtl.CorrelationCalc

/**
 *
 * @author victor
 */
class AssocCorr {
    CorrelationCalc corrCalc
    def scores = [:]
    def sortedItems
    
    static def squareMeasure = {val-> (val==null) ? 0.0 : val*val }
    static def absMeasure = {val-> (val==null) ? 0.0 : Math.abs(val) }
    
    protected def scoreComp = [ compare: {a,b-> scores[b]<=>scores[a] } ] as Comparator
    protected def labelComp = [ compare: {a,b-> Math.abs(corrCalc.corrValues[b])<=>Math.abs(corrCalc.corrValues[a]) } ] as Comparator
    
    static final String OUT_SUFF_SNP = '_assoc_snp.txt'
    static final String OUT_SUFF_ISO = '_assoc_iso.txt'
    static final String WIG_PREF = 'isoform_'
    
    /**
     *
     */
    static def perform(args){
        //read outdir arg
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        def corrFile = outdir+Main.CORREL_FILE
        
        //parse correlation results file
        CorrelationCalc corr = Main.parseCorrFile(corrFile)
       
        //calc assoiation for snps and isoforms
        def assocObj = new AssocCorr(corrCalc:corr)
        
        [OUT_SUFF_SNP, OUT_SUFF_ISO].each{ suff->
            assocObj.calcScores(squareMeasure, suff)
            assocObj.sortScores(suff)

            def outFile = corrFile.substring(0, corrFile.lastIndexOf('.'))
            assocObj.writeResult(outFile+suff, suff)
        }
        
        //generate wig files
        double pvalueThr = 0.001
        double corrThr = 0.4
        generateWigFiles(outdir, pvalueThr, corrThr)
        
    }
    
    /**
     *
     */
    protected def writeResult(file, suff){
        def writer = new PrintWriter(file)
        def xtractAssoc
        
        // write header
        if( suff==OUT_SUFF_SNP ){
            writer.println("snp\tscore\tnum_isos\tsorted_isos")
            xtractAssoc = {lbl-> lbl.substring(lbl.indexOf(':')+1)}
        }
        else{
            writer.println("isoform\tscore\tnum_snps\tsorted_snps")
            xtractAssoc = {lbl-> lbl.substring(0, lbl.indexOf(':'))}
        }
        
        sortedItems.each{ item->
            def list = getSortedAssoc(item, suff)
            //write item, score and num_associated
            writer.print("${item}\t${scores[item]}\t${list.size()}\t")
            //write sorted associated
            list.each{ writer.print("${xtractAssoc(it)}:${Math.abs(corrCalc.corrValues[it])},") }
            writer.print("\n")
        }
        
        writer.close()
    }
    
    /**
     *
     */
    protected def getSortedAssoc(item, suff){
        //builds a list of labels
        def list
        
        if( suff==OUT_SUFF_SNP )
            list = corrCalc.isoforms.collect{corrCalc.buildLabel(item,it)}
        else
            list = corrCalc.snps.collect{corrCalc.buildLabel(it,item)}
            
        list = list.findAll{corrCalc.corrValues[it]!=null} as List
        
        Collections.sort(list, labelComp)
        return list
    }
    
    /**
     *
     */
    protected def sortScores(suff){
        if( suff==OUT_SUFF_SNP )
            sortedItems = corrCalc.snps as List
        else
            sortedItems = corrCalc.isoforms as List 
        
        Collections.sort(sortedItems, scoreComp)
    }
    
    /**
     *
     */
    protected def calcScores(measure, suff){
        if( suff==OUT_SUFF_SNP )
            corrCalc.snps.each{ scores[it] = snpScore(it, measure) }
        else
            corrCalc.isoforms.each{ scores[it] = isoScore(it, measure) }
    }
    
    /**
     *
     */
    protected double snpScore(snp, measure){
        return corrCalc.isoforms.sum{ measure( corrCalc.corrValues[corrCalc.buildLabel(snp,it)] ) }
    }
    
    /**
     *
     */
    protected double isoScore(iso, measure){
        return corrCalc.snps.sum{ measure( corrCalc.corrValues[corrCalc.buildLabel(it,iso)] ) }
    }
    
    /**
     *
     */
    def generateWigFiles(outdir, double pvalueThr, double corrThr){
        def selIsoforms = [] as TreeSet
        
        //select isoforms that will be drawn
        corrCalc.snps.each{ snp->
            corrCalc.isoforms.each{ iso->
                def lbl = corrCalc.buildLabel(snp,iso)
                Double pval = corrCalc.corrPValues[lbl]
                Double corr = corrCalc.corrValues[lbl]
                
                if(  pval!=null && corr!=null && pval<pvalueThr ){
                    if( Math.abs(corr) > corrThr )
                        selIsoforms << iso
                }
            }
        }
        
        //load snps coordinates
        def snps = loadMapFile(outdir+SNPManager.MERGED_PED+'.map')
        
        selIsoforms.each{ iso-> generateWig(iso, snps, outdir) }
    }
    
    /**
     *
     */
    protected def loadMapFile(file){
        def reader = Utils.createReader(file)
        def list = []
        
        reader.splitEachLine("\\s"){ toks->
            list << (new SNPData(chr:toks[0], id:toks[1], position:(toks[3] as Integer) ))
        }
        
        reader.close()
        
        return list
    }
    
    /**
     * generates wig file
     */
    protected def generateWig(iso, snps, outdir){
        def file = outdir+WIG_PREF+iso+'.wig'
        def writer = new PrintWriter(file)
        
        //write header
        writer.println("track type=wiggle_0 name=\"${iso}_corr\" description=\"${iso} correlations\"")
        
        def wigDecl = "variableStep chrom=${snps[0].chr}"
        
        snps.each{ snp->
            def val = corrCalc.corrValues[corrCalc.buildLabel(snp.id,iso)]
            if( val!=null ){
                writer.println(wigDecl)
                writer.println("${snp.position} ${val}")
            }
        }
        
        writer.close()
    }
    
}

