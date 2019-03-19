#!/usr/bin/env groovy

import org.jets3t.service.S3Service
import org.jets3t.service.S3ServiceException
import org.jets3t.service.security.AWSCredentials
import org.jets3t.service.impl.rest.httpclient.RestS3Service
import org.jets3t.service.model.S3Object
import org.jets3t.service.model.S3Bucket
import org.jets3t.service.utils.MultipartUtils

/**
 * Example:
 * groovy -cp libjar/jets3t-0.9.0.jar:libjar/httpclient-4.2-beta1.jar:libjar/httpcore-4.2-beta1.jar:libjar/commons-codec-1.4.jar:libjar/commons-logging-1.1.1.jar:libjar/java-xmlbuilder-0.4.jar geuvadis_download.groovy <group> <output directory> <S3 output folder>
 */

if( !this.args || this.args.length!=3 ){
    println 'Download geuvadis sequences from ENA'
    println 'usage: geuvadis_download <group> <output directory> <S3 output folder>'
    return 1
}


String AWSAccessKeyId = ''
String SecretAccessKey = ''

if( !AWSAccessKeyId || !SecretAccessKey ){
    println 'AWS credentials not setted!'
    return 1
}

long maxSizeForAPartInBytes = 256L * 1024L * 1024L //256 MB
S3Service s3Service

AWSCredentials awsCredentials =  new AWSCredentials(AWSAccessKeyId, SecretAccessKey);
s3Service = new RestS3Service(awsCredentials)

/**
 * Upload a large file to S3 (using multipart API)
 */
def uploadLargeFile = {bucket, file, object ->
    S3Object largeFileObject 

    try{
        largeFileObject = new S3Object(new File(file))
    } catch( FileNotFoundException e){
        return null
    }

    largeFileObject.setName(object)

    MultipartUtils mpUtils = new MultipartUtils(maxSizeForAPartInBytes)

    try{
        mpUtils.uploadObjects( bucket, s3Service, 
            [largeFileObject],
            null // eventListener : Provide one to monitor the upload progress
        )
    } catch(e){
         println "Error uploading to ${bucket+object}"
         return false
    }

    return true
}

// constants
int GROUP_IDX = 6
int FASTQ_IDX = 28
String WGET_COMMAND = "wget -nv --no-proxy"
String infoFileUrl = 'http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt'
String infoFile = infoFileUrl.substring( infoFileUrl.lastIndexOf('/')+1 )

// download closure
def download = {urlSrc, urlDst=null ->
        def target = urlDst ? "-O ${urlDst}" : ''
        def p = "${WGET_COMMAND} ${target} ${urlSrc}".execute()
        
        if(p.waitFor()!=0){
            "rm ${urlDst ?: ''}".execute()//remove empty output file
            return false
        }
        
        return true
    }

// get commandline args
String group = this.args[0].toUpperCase()
String dir = this.args[1]
String s3Out = this.args[2]
String s3Dir = ''
String bucket = s3Out

if( s3Out.indexOf('/')>=0 ){
    bucket = s3Out.substring(0, s3Out.indexOf('/'))
    s3Dir = s3Out.substring(s3Out.indexOf('/')+1)
}

if( !dir.endsWith('/') )
    dir += '/'
    
if( !s3Dir.endsWith('/') )
    s3Dir += '/'

// download geuvadis info file
download(infoFileUrl, dir+infoFile)

int cont = 0
new File(dir+infoFile).splitEachLine("\t"){ toks->
    if( toks[GROUP_IDX]==group ){
        cont++
    }
}

println "Found ${cont} fastq files to download."

//  parse sdrf file and filter by group
new File(dir+infoFile).splitEachLine("\t"){ toks->
    if( toks[GROUP_IDX]==group ){
        //for each subject in group download reads and upload them to S3
        def readFile = toks[FASTQ_IDX].substring( toks[FASTQ_IDX].lastIndexOf('/')+1 )
        
        println "Downloading ${toks[FASTQ_IDX]}, group ${group} [${new Date()}]"
        if( !download(toks[FASTQ_IDX], dir+readFile) )
            println "Error downloading ${toks[FASTQ_IDX]}"
        
        println "Uploading ${dir+readFile}, to ${bucket+s3Dir} [${new Date()}]"
        uploadLargeFile(bucket, dir+readFile, s3Dir+readFile)
             
        //clean up
        "rm -f ${dir+readFile}".execute().waitFor()
    }
}

return 0
