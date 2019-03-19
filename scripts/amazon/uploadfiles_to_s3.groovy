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
 * groovy -cp libjar/jets3t-0.9.0.jar:libjar/httpclient-4.2-beta1.jar:libjar/httpcore-4.2-beta1.jar:libjar/commons-codec-1.4.jar:libjar/commons-logging-1.1.1.jar:libjar/java-xmlbuilder-0.4.jar uploadfiles_to_s3.groovy <list_files> <tmp dir> <bucket>
 */

if( !this.args || this.args.length!=3 ){
    println 'Download files and upload to S3'
    println 'usage: uploadfiles_to_s3 <list_files> <tmp dir> <bucket>'
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

String WGET_COMMAND = "wget -nv --no-proxy"

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
    

String file = this.args[0]
String dir = this.args[1]
String s3Out = this.args[2]
String s3Dir = ''
String bucket = s3Out

if( s3Out.indexOf('/')>=0 ){
    bucket = s3Out.substring(0, s3Out.indexOf('/'))
    s3Dir = s3Out.substring(s3Out.indexOf('/')+1)
}

if( !s3Dir.endsWith('/') )
    s3Dir += '/'
    
if( !dir.endsWith('/') )
    dir += '/'
    
new File(file).eachLine{ line->
    def finalFile = line.substring( line.lastIndexOf('/')+1 )
    
    println "Downloading ${line} [${new Date()}]"
    if( !download(line, dir+finalFile) )
        println "Error downloading ${line}"

    println "Uploading ${dir+finalFile}, to ${bucket+s3Dir} [${new Date()}]"
    uploadLargeFile(bucket, dir+finalFile, s3Dir+finalFile)

    //clean up
    "rm -f ${dir+finalFile}".execute().waitFor()
}

println "Finished ${new Date()}"

return 0
