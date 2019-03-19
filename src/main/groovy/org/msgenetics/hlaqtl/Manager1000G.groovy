/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl


import org.ngsutils.Utils
import org.awsjoblauncher.storage.*

/**
 * 
 * @author victor
 */
class FastqRead {
    String url
    String urlPair
    String sample
    String population
    String comment
    String runId
    String analysisGrp // low coverage, high coverage, exon targetted and exome
    Integer insertSize
    boolean paired
    Long readCount
    Long baseCount
}

/**
 * 
 * 1000genomes manager
 * 
 * @author victor
 */
class Manager1000G {
	
    static final String _1000G_FTP_BASE = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/'
    
    //seq index file fields
    static final int SQI_FASTQ_FILE = 0
    static final int SQI_RUN_ID = 2
    static final int SQI_SAMPLE_NAME = 9 // subject id
    static final int SQI_POPULATION = 10 // group (CEU, FIN, ...)
    static final int SQI_INSERT_SIZE = 17
    static final int SQI_LIBRARY_LAYOUT = 18 // this can be either PAIRED or SINGLE 
    static final int SQI_PAIRED_FASTQ = 19
    static final int SQI_COMMENT = 22
    static final int SQI_READ_COUNT = 23
    static final int SQI_BASE_COUNT = 24
    static final int SQI_ANALYSIS_GROUP = 25 // low coverage, high coverage, exon targetted and exome
    
    static final String SQI_SUPR_COMM = 'SUPPRESSED IN ARCHIVE'
    static final String SQI_FILENAME = 'sequence.index'
    
    /**
     *
     */
    static def prepareReads(args, subjectsMap){
        //read outdir arg
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        String bucket = Main.readOption(args, Main.OPT_S3OUT)
        
        // 1 - download and parse sequence.index
        if( !Utils.download( _1000G_FTP_BASE+SQI_FILENAME, outdir+SQI_FILENAME ) ){
            System.err.println("Error downloading ${SQI_FILENAME}... Aborting")
            System.exit(1)
        }
        
        Manager1000G instance = new Manager1000G()
        
        def readList = instance.parseSeqIndex(outdir+SQI_FILENAME, subjectsMap.keySet(), null)
                
        // 2 - filter and download low coverage reads only
        def filterReads = readList.findAll{ it.analysisGrp.contains('low coverage') }
        
        downloadAndUploadReads(filterReads, outdir, bucket)
        
        // 3 - upload sequence.index to s3
        uploadFiles(outdir, bucket)
    }
    
    /**
     *
     */
    List parseSeqIndex(file, subjectIds, groups) {
        def reader = Utils.createReader(file)
        def fastqReads = []
        
        def formatLong = { str-> try { return str as Long } catch(e){ return null } }
        def formatInt = { str-> try { return str as Integer } catch(e){ return null } }
        
        reader.readLine()//skip header
        reader.splitEachLine("\\t"){ toks->
            boolean keep = false
            
            if( subjectIds && (toks[SQI_SAMPLE_NAME] in subjectIds) )
                keep = true
            else if( groups && (toks[SQI_POPULATION] in groups) )
                keep = true
                
            if( keep && toks[SQI_COMMENT].contains(SQI_SUPR_COMM) )
                keep = false //skip suppressed files
                
            if( keep ){
                fastqReads << new FastqRead (
                    url: _1000G_FTP_BASE+toks[SQI_FASTQ_FILE],
                    urlPair: toks[SQI_PAIRED_FASTQ] ? _1000G_FTP_BASE+toks[SQI_PAIRED_FASTQ] : '',
                    sample: toks[SQI_SAMPLE_NAME],
                    population: toks[SQI_POPULATION],
                    comment: toks[SQI_COMMENT],
                    runId: toks[SQI_RUN_ID],
                    insertSize: formatInt(toks[SQI_INSERT_SIZE]),
                    paired: (toks[SQI_LIBRARY_LAYOUT]=='PAIRED'),
                    readCount: formatLong(toks[SQI_READ_COUNT]),
                    baseCount: formatLong(toks[SQI_BASE_COUNT]),
                    analysisGrp: toks[SQI_ANALYSIS_GROUP] )
            }
        }
        
        reader.close()
        return fastqReads
    }
    
    /**
     *
     */
    def downloadAndUploadReads(reads, outdir, bucket) {
        String s3Dir = ''
        String s3Out = bucket

        if( s3Out.indexOf('/')>=0 ){
            bucket = s3Out.substring(0, s3Out.indexOf('/'))
            s3Dir = s3Out.substring(s3Out.indexOf('/')+1)
        }

        if( !s3Dir.endsWith('/') )
            s3Dir += '/'
            
        reads.each{ read->
            String name = read.url.substring( read.url.lastIndexOf('/')+1 )
            
            if( !Utils.download( read.url, outdir+name ) )
                println "Error downloading ${read.url}"
                
            if( !StorageManager.uploadToS3(bucket, outdir+name, s3Dir+name) )
                println "Error uploading to ${s3Dir+name}"
                    
            //clean up
            "rm -f ${outdir+name}".execute().waitFor()
        }
    }
   
    /**
     *
     */
    def uploadFiles(dir, bucket) {
        String s3Dir = ''
        String s3Out = bucket

        if( s3Out.indexOf('/')>=0 ){
            bucket = s3Out.substring(0, s3Out.indexOf('/'))
            s3Dir = s3Out.substring(s3Out.indexOf('/')+1)
        }

        if( !s3Dir.endsWith('/') )
            s3Dir += '/'
    
        new File(dir).eachFile{ file->
            if( file.isFile() ){
                if( !StorageManager.uploadToS3(bucket, dir+file.name, s3Dir+file.name) )
                    println "Error uploading to ${s3Dir+file.name}"
            }
        }
    }
    
}

