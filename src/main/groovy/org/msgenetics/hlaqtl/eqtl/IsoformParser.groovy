/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.awsjoblauncher.storage.*
import org.awsjoblauncher.postprocess.CommConstants
import org.ngsutils.annotation.*
import org.ngsutils.Utils
import org.ngsutils.diff.IsoformExpData
import org.ngsutils.annotation.biojava.GFFWriterLight

/**
 * Class that represents the Isoform expression data of an individual
 * 
 * @author victor
 */
class IsoformParser {
    def subject // a Subject object
    String bucket
    String workDir
    String cuffDirPrefix
    
    def isoforms = [:] //map with key:<Isoform id> and value:IsoformExpData
    
    static final String ISO_GTF_FILE = CommConstants.GTF_OUTPUT+'.gz'
        
    static final String CUFF_ISO_PRE = 'CUFF.'
    static final String NOVEL_ISO_PRE = 'NEW.'
    
    static final int TRK_ISO_ID = 2 //isoform id field in subject isoforms tracking file
    static final int CUFF_ISO_ID = 0//isoform id field in cufflinks tracking file
    
    /**
     *
     */
    protected def parseRegion(String locusStr, String exprdir, selIso=null) {
        if( locusStr==null ) {
            // keep all isoforms (trans-eqtl studies)
            return parseRegion(null, 0, 0, exprdir, selIso)
        }
        
        def locus = (locusStr =~ Utils.locusRegex)
        parseRegion(locus[0][1], locus[0][2] as Integer, locus[0][3] as Integer, exprdir, selIso)
    }
    
    /**
     *
     */
    protected String outTrkFile(dir) {
        "${dir}${subject.id}_${IsoformExpData.ISO_TRACK_FILE}"
    }
    
    /**
     *
     */
    protected String outGtfFile(dir) {
        "${dir}${subject.id}_${ISO_GTF_FILE}"
    }
    
    /**
     *
     */
    protected def downloadExprData() {
        def cuffDir = cuffDirPrefix+subject.id+'/'
        def trackingFile = cuffDir+IsoformExpData.ISO_TRACK_FILE
        def gtfFile = cuffDir+ISO_GTF_FILE
        
        try {
            if( (new File(outTrkFile(workDir))).exists() ) {
                println "Subject [${subject.id}]: skip downloading ${trackingFile}"
            }
            else if( S3Manager.downloadFile( bucket, trackingFile, outTrkFile(workDir))==null ) {
                ///throw new StorageException("Error downloading ${trackingFile}")
                println "Subject [${subject.id}]: Error downloading ${trackingFile}"
                return null
            }
            
            if( (new File(outGtfFile(workDir))).exists() ) {
                println "Subject [${subject.id}]: skip downloading ${gtfFile}"
            }
            else if( S3Manager.downloadFile( bucket, gtfFile, outGtfFile(workDir))==null ){
                ///throw new StorageException("Error downloading ${gtfFile}")
                println "Subject [${subject.id}]: Error downloading ${gtfFile}"
                return null
            }
        }
        catch(e) {
           println "Subject [${subject.id}]: Exception while downloading files"
           [outTrkFile(workDir),outGtfFile(workDir)].each{ "rm -f ${it}".execute() }
           return null 
        }
        
        return true
    }
    
    /**
     * downloads cufflinks results from bucket/cuffDirPrefix${subject} and parses 
     * gtf and tracking files
     * 
     * @param chr : chromosome to parse or null to keep all isoforms
     */
    protected def parseRegion(String chr, int start, int end, String exprdir, selIso=null) {
        def dir = exprdir
        
        if( !exprdir ){
            //download cufflinks_ref result
            if( !downloadExprData() )
                return null
                
            dir = workDir
        }
        
        if( !(new File(outTrkFile(dir))).exists() ){
            println "Subject [${subject.id}]: Warning ${outTrkFile(dir)} not exists"
            return null
        }
            
        //2 - parse gtf file
        def geneIndex = ((chr!=null) ? new GeneIndex(new File(outGtfFile(dir)), chr) : null)
        
        //3 - parse tracking file
        def reader = Utils.createReader(outTrkFile(dir))
        reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            if( toks[IsoformExpData.STATUS_IDX]=='OK' &&  
                (chr==null || IsoformExpData.isInRegion(chr, start, end, toks[IsoformExpData.LOCUS_IDX]))
                ) {
                if( !selIso || (selIso && (toks[CUFF_ISO_ID] in selIso)) ){
                    //keep isoforms located inside the desired region
                    def isoform = IsoformExpData.parse(toks)
                    isoform.isoform = ((geneIndex!=null) ? geneIndex.getIsoform(isoform.id) : null)
                    isoforms[isoform.id] = isoform
                }
            }
        }
        
        reader.close()
        
        //clean tmp files
        if( !exprdir )
            [outTrkFile(dir), outGtfFile(dir)].each{ "rm -f ${it}".execute().waitFor() }
    }
     
    /**
     *
     * @param list : list of IsoformParser objects, one per subject
     */
    protected static mergeIsoforms(list){
        def mergedIso = [:]//a map of isoforms
        
        def isInMerged = { iso-> mergedIso.find{ it.value.isEqualIso(iso) } }
        
        list.each{ isoParser->
            isoParser.isoforms.each{ id, isoData->
                if( id.startsWith(CUFF_ISO_PRE) ){
                    def merged = isInMerged(isoData.isoform)
                    
                    if( merged ){
                        isoData.mergedId = merged.key
                    }
                    else{
                        isoData.mergedId = "${NOVEL_ISO_PRE}${mergedIso.size()+1}"
                        mergedIso[isoData.mergedId] = isoData.isoform
                    }
                    
                    isoData.isoform.transcriptId = isoData.mergedId
                }
            }
        }
        //reindex isoforms: new key = merged_id
        list.each{ isoParser->
            def merged = isoParser.isoforms.findAll{it.key.startsWith(CUFF_ISO_PRE)}
            
            merged.each{ id, isoData->
                isoParser.isoforms[isoData.mergedId] = isoData
                isoParser.isoforms.remove(id)
            }
        }
    }
    
    /**
     * writes a gff file with the set of merged isoforms
     * 
     * @param list : list of IsoformParser objects, one per subject, previously
     * merged
     */
    static writeIsoforms(list, file){
        def mergedIso = [:] as TreeMap//a map of isoforms
        
        //populate isoform map
        list.each{ isoParser->
            isoParser.isoforms.each{ id, isoData->
                def key = isoData.mergedId ?: isoData.id
                
                if( !mergedIso.containsKey(key) )
                    mergedIso[key] = isoData
            }
        }
        
        //sort isoforms by locus (chr:start-end)
        def sortedList = mergedIso.values().sort(IsoformExpData.locusComparator)
        
        //write isoforms to gff file
        def writer = new PrintWriter(file)
        def gffWriter = new GFFWriterLight(writer)
        
        sortedList.each{ it.isoform?.record(gffWriter) }
        
        writer.close()
    }
    
    
    /**
     * writes readed isoforms (id, merged id etc.) of each subject
     * 
     * @param list : list of IsoformParser objects, one per subject, previously
     * merged
     */
    static writeSubjectIsoforms(list, file){
        def list2 = list.sort{ it.subject.id }
        
        def writer = new PrintWriter(file)
        
        list2.each{ isoParser->
            def subjId = isoParser.subject.id
            
            isoParser.isoforms.values().sort(IsoformExpData.locusComparator).each{
                writer.println(
                "${subjId}\t${it.id}\t${it.mergedId ?: it.id}\t${it.chr}:${it.start}-${it.end}\t${it.fpkm}"
                )
            }
        }
        
        writer.close()
        
        "gzip -f ${file}".execute().waitFor()
    }
    
    
    /**
     * loads and prepares expression data of each subject
     * 
     * @param subjects : a list of Subject object
     * @param exprdir : directory where to find expression data
     * @param selIso : selected isoforms Ids set
     * @return a list of IsoformParser objects (one per subject)
     */
     static loadExpressionData(subjects, resBucket, cuffDirPrefix, outdir, locus, exprdir, selIso = null){
        def listIsoData = subjects.collect{ 
            new IsoformParser( subject:it, bucket:resBucket, workDir:outdir, 
                cuffDirPrefix:cuffDirPrefix ) 
        }
        
        listIsoData.each{ it.parseRegion(locus, exprdir, selIso) }
        IsoformParser.mergeIsoforms(listIsoData)
        
        return listIsoData
    }
    
    /**
     *
     */
    static downloadExpressionData(subjects, resBucket, cuffDirPrefix, outdir){
        def listIsoData = subjects.collect{ 
            new IsoformParser( subject:it, bucket:resBucket, workDir:outdir, 
                cuffDirPrefix:cuffDirPrefix ) 
        }
        
        listIsoData.each{ it.downloadExprData() }
        
        return true
    }
    
    /**
     *
     */
    static protected IsoformExpData newIsoformExpData(toks) {
        def locus = (toks[3] =~ Utils.locusRegex)
        
        return new IsoformExpData( id:toks[TRK_ISO_ID], mergedId:toks[TRK_ISO_ID],
                chr:locus[0][1], start:(locus[0][2] as Integer), 
                end:(locus[0][3] as Integer), fpkm:(toks[4] as Float) )
    }
    
    /**
     * parses isoforms of each subject
     * 
     * @param file : iso_merged.tracking file
     * @param subjsMap  : map with key=subjectId and value=Main.Subject object
     * @return a list of IsoformParser objects (one per subject)
     */
    static parseSubjectIsoforms(file, subjsMap, selIso=null) {
        def listIsoData = []
        def currSubj = null
        def currIsoData = null
        
        def reader = Utils.createReader(file)
        
        reader.splitEachLine("\\s"){ toks->
            //subj iso iso_merg locus fpkm
            if( toks[0]!=currSubj ){
                currSubj = toks[0]
                currIsoData = new IsoformParser(subject:subjsMap[currSubj])
                listIsoData << currIsoData
            }
            
            if( !selIso || (selIso && (toks[TRK_ISO_ID] in selIso)) ){
                currIsoData.isoforms[toks[TRK_ISO_ID]] = newIsoformExpData(toks)
            }
        }
        
        reader.close()
        
        return listIsoData
    }
    
    /**
     * parses distinct isoforms in file 
     * 
     * @return a map of IsoformExpData objects (key=iso_id)
     */
    static parseIsoforms(file, map = [:]) {
        def reader = Utils.createReader(file)
        
        reader.splitEachLine("\\s"){ toks->
            if( !map.containsKey(toks[TRK_ISO_ID]) ) {
                map[toks[TRK_ISO_ID]] = newIsoformExpData(toks)
            }
        }
        
        reader.close()
        
        return map
    }
    
    /**
     * parse a isoforms filtered file
     * 
     * @return a Set of isoforms IDs
     */
    static parseIsoFilterFile(file){
        def isoforms = [] as TreeSet
        
        // parse isoforms filtered file
        def reader = new File(file).newReader()
        reader.readLine()//skip header
        reader.splitEachLine("\t"){ toks-> isoforms << toks[0] }
        reader.close()
        
        return isoforms
    }
    
    /**
     * parse a isoforms filtered file
     * 
     * @return a Map with pairs key=isoformId and value=gene_name
     */
    static parseIsoFilterFileGenes(file){
        def isoforms = [] as TreeMap
        
        // parse isoforms filtered file
        def reader = new File(file).newReader()
        reader.readLine()//skip header
        reader.splitEachLine("\t"){ toks-> isoforms[toks[0]] = toks[6] }
        reader.close()
        
        return isoforms
    }
    
}

