#!/usr/bin/env groovy

/**
 * generates commands for eqtl trans calculation in separate files
 */

int THREADS = 6
int SNP_SPLIT_NUM = 300
String SNP_SPLIT_PREF = 'list_snps_part.'
String EQTL_COMM_PREF = 'eqtl_commands.'
boolean QSUB = true
String email = 'vpoten@gmail.com'
int memory = 62

// directories and options
String javaBin = 'java'
String jarFile = '/home/users/ipb/lindo/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar'
String javaOpts = "-Xmx${memory}g"
String awsProps = '/home/users/ipb/lindo/HlaQtl/ngsengine_eqtltrans.properties'
String outDir = '/home/users/ipb/lindo/eqtl_missense_out/'
String exprDataDir = '/home/users/ipb/lindo/cufflinks_ensOnly/'
String groups = 'CEU,GBR,FIN,TSI'
String cuffPre = 'cufflinks_ensOnly_'
String vcfDir = '/home/users/ipb/lindo/vcfFiles/'
String snpsFile = '/home/users/ipb/lindo/missense_snps.freq_filter'

// split snpsFile
"split -l ${SNP_SPLIT_NUM} ${snpsFile} ${outDir+SNP_SPLIT_PREF}".execute().waitFor()

// get splitted files
def snpFiles = ((new File(outDir)).listFiles() as List).findAll{ it.name.contains(SNP_SPLIT_PREF) }.sort{it.name}

def commands = []
(1..THREADS).each{ commands << [] }

snpFiles.eachWithIndex{ file, i->
    int idx = i%THREADS
    def dir = 'eqtl_trans_'+file.name.substring( file.name.lastIndexOf('.')+1 )
    
    //create dir
    "mkdir ${outDir+dir}".execute().waitFor()
    
    def comm = "${javaBin} ${javaOpts} -jar ${jarFile} eqtl --props=${awsProps} " +
        "--output=${outDir+dir} --group=${groups} --cuffPre=${cuffPre} " +
        "--exprData=${exprDataDir} --vcf=${vcfDir} --transEqtlSnps=${file.absolutePath} "+
        "> ${outDir+dir}.out"
    
    commands[idx] << comm
}

//write commands to files
commands.eachWithIndex{ list, i->
    def writer = new PrintWriter("${outDir}${EQTL_COMM_PREF}${i}")
    
    if( QSUB ){
        //write qsub header
        writer.println('#!/bin/bash')
        writer.println("#PBS -N eqtlTrans_${i}")
        writer.println("#PBS -M ${email}")
        writer.println('#PBS -m abe')
        writer.println("#PBS -l nodes=1:ppn=6,mem=${memory}gb")
    }
    
    list.each{ writer.println(it) }
    writer.close()
}

return 0