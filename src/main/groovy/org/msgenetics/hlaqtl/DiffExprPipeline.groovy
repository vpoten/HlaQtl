/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.awsjoblauncher.storage.StorageManager
import org.awsjoblauncher.postprocess.CommConstants

class DESubject {
    String id
    String group
    String condition
    def fastq1
    def fastq2
}

/**
 *
 * @author victor
 */
class DiffExprPipeline {
	
    /**
     *
     */
    static perform(args) {
        String expFile = Main.readOption(args, Main.OPT_SUBJECTS)
        int start = Main.getStartPhase(args)
        String value = Main.readOption(args, Main.OPT_VALUE)
        
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        def dirs = Main.checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        //// Options
        ///String bwtIndex = Main.config.org.webngs.bowtieIdx.dir+'hg19'
        String bwtIndex = Main.config.org.webngs.bowtieIdx.dir + 
            'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
        ///String reference = StorageManager.ANNOTATIONS_PRE+'9606/ensGene.txt.gtf.gz'
        String reference = StorageManager.ANNOTATIONS_PRE+'9606/ensembl_hg38_filtered.gff3'
        String topHatRef = null
        ////
        
        println "experiments file: ${expFile}"
        println "output: ${outdir}"
        println "bowtie index: ${bwtIndex}"
        println "reference annotation: ${reference}"
        println "start phase: ${start}"
        
        // read groups and subjects
        def groups = readGroups(expFile)
        
        def jobList = []
        boolean execLocal = false
        
        
        
        if( start==0 ) {
            // FastQC
            groups.values().each{ list-> 
                list.each{ jobList << JobGenerator.fastQC(it.id, StorageManager.UPLOADS_PRE, it.fastq1+it.fastq2) }
            }
        }
        else if( start==1 ) {
            // FastxToolkit
            execLocal = true
            String inputUrl = StorageManager.UPLOADS_PRE
        
            println "Filtering of fastq data in ${inputUrl}"
            def optMap = [minLen:'50', maxLen:'no', trimFirst:'14', qualFilt:'20', percQual:'66']
            println "Filtering parameters:\n${optMap}"
        
            groups.values().each{ list-> 
                list.each{ jobList << JobGenerator.fastxToolkit(it.id, inputUrl, optMap, it.fastq1, it.fastq2) }
            }
        }
        else if( start==2 ){
            // tophat + cufflinks
            groups.values().each{ list->
                list.each{
                    def input = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTXTOOL}_${it.id}/"
                    jobList << JobGenerator.tophat(it.id, input, bwtIndex, topHatRef)
                    
                    input = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${it.id}/"
                    jobList << JobGenerator.cufflinks(it.id, input, reference, true)
                }
            }
        }
        else if( start==3 ){
            // cuffdiff
            groups.each{ group, list->
                def subjIds = list.collect{it.id}
                def jobId = value ? "${value}-g${group}": Integer.toHexString(subjIds.sum().hashCode()).toUpperCase()
                jobList << JobGenerator.cuffdiff(jobId, subjIds, reference)
            }
        }
        
        if( execLocal ){
            Main.execJobs(jobList)
            return
        }
        
        def writer = new PrintWriter( outdir+'diffexpr_commands.sh.txt')
        writeHeader(writer)
        
        jobList.each{ job->
            writer.println("#Job: ${job.jobId}")
            writer.println("mkdir ${StorageManager.replacePrefixs(job.url)}")
            writer.println(StorageManager.replacePrefixs(job.command))
        }
        
        writer.close()
    }
    
    /**
     * print qsub header
     */ 
    static def writeHeader(writer) {
        writer.println("#!/bin/bash")
        writer.println("#")
        writer.println("#PBS -N diffExprPipe_")
        writer.println("#PBS -M <email>")
        writer.println("#PBS -m abe")
        writer.println("#PBS -l nodes=1:ppn=7,mem=8gb")
        writer.println('#')    
    }   
    
    /**
     *
     */
    static readGroups(expFile) {
        def groups = [:]
        
        def reader = new File(expFile).newReader()
        reader.readLine()//skip header
        
        reader.splitEachLine("\t"){ toks->
            //file fields: subject group label fastq_1 fastq_2
            def subj = new DESubject( id:toks[0], group:toks[1], condition:toks[2],
                fastq1:(toks[3].split(',') as List), fastq2:(toks[4].split(',') as List) )
            
            def grp = groups[subj.group]
            
            if( grp==null ) {
                grp = []
                groups[subj.group] = grp
            }
            
            grp << subj
        }
        
        reader.close()
        
        return groups
    }
    
}

