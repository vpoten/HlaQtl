/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.awsjoblauncher.storage.*
import org.awsjoblauncher.postprocess.CommConstants

/**
 *
 * @author victor
 */
class DownJobData {
    
    /**
     * download to local tophat + fastxtoolkit data of each subject
     */
    static def downloadJobData(args, subjects) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
            
        def folders = []
        
        //tophat + fastxtoolkit data
        subjects.each{
            folders  << "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTXTOOL}_${it.key}"
            folders  << "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${it.key}"
        }
        
        folders.each{ url->
            def bucket = StorageManager.getBucket(url)
            def name = StorageManager.removePrefix(url)
            def result = S3Manager.listObjects(bucket, name);
            
            if( result )
                "mkdir ${outdir+name}".execute().waitFor()
            else
                println "No files inside ${url} folder; bucket [${bucket}]"
            
            result.each{
                println "Downloading ${name+'/'+it}; bucket [${bucket}]"
                try {
                    if( !StorageManager.downloadFromS3(bucket, name+'/'+it, outdir+name+'/'+it) )
                        println "Error downloading ${name+'/'+it}"
                }
                catch(e){
                    println "Error downloading ${name+'/'+it} (due an Exception)"
                }
            }
        }
    }
	
    /**
     * download to local eqtl result in output bucket
     */
    static def downloadEqtlResult(args) {
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        def locus = Main.readOption(args, Main.OPT_LOCUS)
        def groups = Main.getOptGroups(args)
        def cuffDirPrefix = Main.readOption(args, Main.OPT_CUFFPRE)
        
        //print readed options
        println "${Main.LINE_PRE}Options present:"
        println "locus: ${locus}"
        println "groups: ${groups}"
        println "output: ${outdir}"
        println "cuffdiff pref: ${cuffDirPrefix}"
        
        String objName = Main.compEqtlResObject(locus, cuffDirPrefix, groups)
        String url = StorageManager.OUTPUT_PRE+objName
        def bucket = StorageManager.getBucket(url)
        
        def files = StorageManager.listFiles(url)
        
        files = files?.findAll{ it.endsWith('.tar.gz') }
        
        if( !files ){
            println "No *.tar.gz files inside ${url}"
            return
        }
        
        try {
            StorageManager.createOutDir(outdir+objName)
        } catch(e) {
            println "Error: cannot create ${outdir+objName} folder"
            return
        }
        
        files.each{
            println "Downloading ${it} to ${outdir+objName}"
            
            if( !StorageManager.downloadFromS3(bucket, objName+'/'+it, outdir+objName+'/'+it) )
                println "Error downloading ${it}"
        }
    }
    
}

