/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.ngsutils.Utils
import org.msgenetics.hlaqtl.eqtl.*
import org.ngsutils.eqtl.CorrelationCalc

/**
 *
 * @author victor
 */
class DebugEqtl {
    
    public static final String SUBJ_DEBUG_FILE = 'debug_subjects.txt'
    public static final String NEWISO_DEBUG_FILE = 'debug_newisos.txt'
    public static final String NEWISOCORR_DEBUG_FILE = 'debug_newisos_corr.txt'
    
    /**
     *
     */
    static def perform(args, subjects){
        //read outdir arg
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        //inspect isoforms tracking
        def listIsoData = IsoformParser.parseSubjectIsoforms(outdir+Main.MERG_TRACK_FILE, subjects)
        
        //write subjects debug
        def writer = new PrintWriter(outdir+SUBJ_DEBUG_FILE)
        
        listIsoData.each{
            def newIsos = it.isoforms.findAll{ it.key.startsWith(IsoformParser.NOVEL_ISO_PRE) }
            
            writer.print("${it.subject.id}\t${newIsos.size()}\t")
            
            newIsos.each{ id, iso-> writer.print("${id}=${iso.fpkm},") }
            
            writer.print("\n")
        }
         
        writer.close()
        
        //write new isos debug
        writer = new PrintWriter(outdir+NEWISO_DEBUG_FILE)
        
        //gather new isoforms
        def isoMap = [:] as TreeMap
        
        listIsoData.each{ isoData->
            isoData.isoforms.each{key, val-> 
                if( val.mergedId.startsWith(IsoformParser.NOVEL_ISO_PRE) ){
                    if( !isoMap.containsKey(val.mergedId) )
                        isoMap[val.mergedId] = []
                        
                    isoMap[val.mergedId] << isoData
                }
            }
        }
        
        isoMap.each{ iso, listData->
            writer.print("${iso}\t${listData.size()}\t")
            
            listData.each{ writer.print("${it.subject.id}=${it.isoforms[iso].fpkm},") }
            
            writer.print("\n")
        }
        
        writer.close()
        
        //parse correlation results file
        CorrelationCalc corrCalc = Main.parseCorrFile(outdir+Main.CORREL_FILE)
        
        //print corr result of new isoforms
        writer = new PrintWriter(outdir+NEWISOCORR_DEBUG_FILE)
        double pvalue = 0.05
        writer.println("isoform\tnum_subjs\tsnp\tcorrelation\tpvalue")
        
        isoMap.each{ iso, listData->
            //get corr snps-iso with pvalue>threshold
            def snps = corrCalc.snps.findAll{ 
                def lbl = corrCalc.buildLabel(it,iso)
                def val = corrCalc.corrPValues[lbl]
                (val!=null && val<pvalue )
            }
            //print iso snp corr pvalue
            snps.each{ snp->
                def lbl = corrCalc.buildLabel(snp,iso)
                writer.println("${iso}\t${listData.size()}\t${snp}\t${corrCalc.corrValues[lbl]}\t${corrCalc.corrPValues[lbl]}")
            }
        }
        
        writer.close()
    }
    
    /**
     *
     */
    protected static getSnpsCorr(CorrelationCalc corrCalc, iso, double pvalue){
        
    }
    
}

