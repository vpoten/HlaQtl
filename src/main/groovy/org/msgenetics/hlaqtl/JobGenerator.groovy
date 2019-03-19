/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.msgenetics.hlaqtl.Subject
import org.awsjoblauncher.postprocess.*
import org.awsjoblauncher.storage.StorageManager
import org.ngsutils.Utils

/**
 *
 * @author victor
 */
class JobGenerator {
    
    static int nThreads = 6
    static int pairDist = 50
    static String organismId = '9606'
    
    //constants
    static def COMMNAME_PHLAT = 'phlat'
    static def JOB_TYPE_PHLAT = 'PHLAT'
    static def PHLAT_FOLDER = "phlat-release/";
    
    /**
     * 
     * @param inputUrl: url base where to find CEU fastq
     */
    static JobDataBean fastQC(Subject subject, inputUrl){
        def inputs = ['1','2'].sum{ "${inputUrl}${subject.lanes[0]}_${it}.fastq.gz " }
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTQC}_${subject.id}"
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_FASTQC}_${subject.id}",
            jobType: CommConstants.JOB_TYPE_FASTQC, 
            command: PostFastQC.generateCommand(outdir, inputs), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    static JobDataBean fastQC(String subjId, inputUrl, fastqList){
        def inputs = fastqList.sum{ "${inputUrl}${it} " }
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTQC}_${subjId}"
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_FASTQC}_${subjId}",
            jobType: CommConstants.JOB_TYPE_FASTQC, 
            command: PostFastQC.generateCommand(outdir, inputs), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    
    /**
     * 
     * @param inputUrl: url base where to find CEU fastq
     * @param optMap : {minLen, maxLen, trimFirst, qualFilt, percQual}
     */
    static JobDataBean fastxToolkit(Subject subject, inputUrl, optMap){
        def inputs = ['1','2'].sum{ "-${it} ${inputUrl}${subject.lanes[0]}_${it}.fastq.gz " }
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTXTOOL}_${subject.id}"
        def options = PostFastxTool.generateOpts( optMap['minLen'], 
            optMap['maxLen'], optMap['trimFirst'], 
            optMap['qualFilt'], optMap['percQual'] )
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_FASTXTOOL}_${subject.id}",
            jobType: CommConstants.JOB_TYPE_FASTXTOOL, 
            command: PostFastxTool.generateCommand(inputs, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    static JobDataBean fastxToolkit(String subjId, inputUrl, optMap, fastq1, fastq2){
        def inputs = "-1 ${fastq1.sum{",${inputUrl}${it}"}.substring(1)} -2 ${fastq2.sum{",${inputUrl}${it}"}.substring(1)}"
        
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTXTOOL}_${subjId}"
        def options = PostFastxTool.generateOpts( optMap['minLen'], 
            optMap['maxLen'], optMap['trimFirst'], 
            optMap['qualFilt'], optMap['percQual'] )
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_FASTXTOOL}_${subjId}",
            jobType: CommConstants.JOB_TYPE_FASTXTOOL, 
            command: PostFastxTool.generateCommand(inputs, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    
    /**
     *
     * @param inputUrl: url base where to find CEU fastq or fastx-toolkit result url
     * @param bwtIndex: bowtie index base
     */
    static JobDataBean tophat(Subject subject, inputUrl, bwtIndex){
        def inputs = ['1','2'].sum{ "${inputUrl}${subject.lanes[0]}_${it}.fastq.gz " }
        
        if( inputUrl.contains(CommConstants.COMMNAME_FASTXTOOL+'_') )//if the input is fastx-toolkit
            inputs = ['1','2'].sum{ "${inputUrl}${CommConstants.FASTXTOOL_OUT}${it}.fastq.gz " }
        
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${subject.id}"
        def options = "-p ${nThreads} -r ${pairDist}"
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_TOPHAT}_${subject.id}",
            jobType: CommConstants.JOB_TYPE_TOPHAT, 
            command: PostTophat.generateCommand(outdir, inputs, bwtIndex, options), 
            url: outdir, genome: organismId, genOutput: true )
        
        return job
    }
    
    static JobDataBean tophat(String subjId, inputUrl, bwtIndex, reference){
        def inputs = ['1','2'].sum{ "${inputUrl}${CommConstants.FASTXTOOL_OUT}${it}.fastq.gz " }
        
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${subjId}"
        def options = "-p ${nThreads} -r ${pairDist}"
        
        if( reference ){
            options += " -G ${reference} --no-novel-juncs"
        }
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_TOPHAT}_${subjId}",
            jobType: CommConstants.JOB_TYPE_TOPHAT, 
            command: PostTophat.generateCommand(outdir, inputs, bwtIndex, options), 
            url: outdir, genome: organismId, genOutput: true )
        
        return job
    }
    
    
    /**
     * 
     * @param inputUrl: tophat result url
     */
    static JobDataBean cufflinks(Subject subject, String inputUrl, String reference, boolean refOnly = false){
        def refPart = ''//reference used (part of output dir name)
        String refPre = ''
        
        if( reference )
            refPre = reference.substring( reference.lastIndexOf('/')+1, reference.indexOf('Gene.') )
        
        if( refOnly )
            refPart ="_${refPre}Only"
        else
            refPart = (reference) ? "_${refPre}" : ''
            
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFLINKS}${refPart}_${subject.id}"
        def inputs = inputUrl+CommConstants.BAM_OUTPUT
        def options = "-p ${nThreads}"
        
        if( refOnly )
            options +=" -G ${reference}"
        else
            options += (reference) ? " -g ${reference}" : ''
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_CUFFLINKS}_${subject.id}",
            jobType: CommConstants.JOB_TYPE_CUFFLINKS, 
            command: PostCufflinks.generateCommand(outdir, inputs, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    static JobDataBean cufflinks(String subjId, String inputUrl, String reference, boolean refOnly = false){
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFLINKS}_${subjId}"
        def inputs = inputUrl+CommConstants.BAM_OUTPUT
        def options = "-p ${nThreads} --no-update-check"
        
        if( refOnly )
            options +=" -G ${reference}"
        else
            options += (reference) ? " -g ${reference}" : ''
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_CUFFLINKS}_${subjId}",
            jobType: CommConstants.JOB_TYPE_CUFFLINKS, 
            command: PostCufflinks.generateCommand(outdir, inputs, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    /**
     * 
     * 
     */
    static JobDataBean cuffdiff(String jobId, List subjIds, String reference) {
        def inputs = subjIds.sum{
            "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${it}/${CommConstants.BAM_OUTPUT} " 
        }
        
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFDIFF}_${jobId}"
        def options = "-p ${nThreads} --no-update-check -L "+subjIds.sum{','+it}.substring(1)
        options += " --dispersion-method blind"//no replicates
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_CUFFDIFF}_${jobId}",
            jobType: CommConstants.JOB_TYPE_CUFFDIFF, 
            command: PostCuffdiff.generateCommand(outdir, inputs, reference, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    
    /**
     * 
     * @param inputList: a list of cufflinks result urls
     * @reference : a url to a gtf reference annotation
     */
    static JobDataBean cuffcompare(Subject subject, List inputList, String reference){
        def outdir = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFCOMPARE}_${subject.id}"
        def outfile = "${outdir}/${CommConstants.COMMNAME_CUFFCOMPARE}_${subject.id}"
        def options = (reference) ? "-r ${reference}" : ''
        
        def inputs = inputList.sum{ it+CommConstants.GTF_OUTPUT+'.gz ' }
        
        def job = new JobDataBean( jobId: "${CommConstants.COMMNAME_CUFFCOMPARE}_${subject.id}",
            jobType: CommConstants.JOB_TYPE_CUFFCOMPARE, 
            command: PostCuffcompare.generateCommand(outfile, inputs, options), 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
    /**
     *
     */
    static String tophatHLA(Subject subject, fastxBase, workDir, boolean createOut=false) {
        def reads = ''
        [1,2].each{ 
            reads += \
        " ${fastxBase}${CommConstants.COMMNAME_FASTXTOOL}_${subject.id}/${CommConstants.FASTXTOOL_OUT}${it}.fastq.gz"
        }
        
        def bwtIndex = "${workDir}${subject.id}/${subject.id}"
        
        def outDir = "${workDir}${CommConstants.COMMNAME_TOPHAT}_${subject.id}"
        
        def options = "-p ${nThreads} -r ${pairDist}"
        
        def gff = "${workDir}${subject.id}/${subject.id}.gff3"
        def extOptions = "--b2-very-sensitive --b2-score-min L,-12,0.0 --b2-np 3"
        extOptions += " -G ${gff} --no-novel-juncs"
        
        if( createOut ) {
            Utils.createDir(outDir)
        }
        
        "${CommConstants.COMMNAME_TOPHAT} ${options} ${extOptions} -o ${outDir} ${bwtIndex} ${reads}"
    }
    
    /**
     *
     */
    static String cufflinksHLA(Subject subject, workDir, boolean createOut=false) {
        def options = "-p ${nThreads}"
        options +=" -G ${workDir}${subject.id}/${subject.id}.gff3"
        
        def outDir = "${workDir}${CommConstants.COMMNAME_CUFFLINKS}_${subject.id}"
        
        def bamInput = "${workDir}${CommConstants.COMMNAME_TOPHAT}_${subject.id}/${CommConstants.BAM_OUTPUT}"
        
        if( createOut ) {
            Utils.createDir(outDir)
        }
        
        "${CommConstants.COMMNAME_CUFFLINKS} ${options} -o ${outDir} ${bamInput}"
    }
    
    /**
     * HLA typing using PHLAT algorythm
     */
    static JobDataBean phlatHLATyping(String subjId, inputUrl, bwtUrl, scriptDir) {
        def inputs = ['1','2'].sum{
        "-${it} ${inputUrl}${CommConstants.COMMNAME_FASTXTOOL}_${subjId}/${CommConstants.FASTXTOOL_OUT}${it}.fastq.gz "
        }
        
        def outdir = "${StorageManager.OUTPUT_PRE}${COMMNAME_PHLAT}_${subjId}"
        def options = "-p ${nThreads} -orientation \"--fr\""
        def phlatdir = scriptDir+PHLAT_FOLDER
        
       
        def phlatCmd = "python2.7 -O ${phlatdir}dist/PHLAT.py ${inputs} -index ${phlatdir}b2folder " +
            "-b2url ${bwtUrl} ${options} -tag ${subjId} -e ${phlatdir} -o ${outdir}"
       
        def job = new JobDataBean( jobId: "${COMMNAME_PHLAT}_${subjId}",
            jobType: JOB_TYPE_PHLAT, 
            command: phlatCmd, 
            url: outdir, genome: organismId, genOutput: false )
        
        return job
    }
    
}

