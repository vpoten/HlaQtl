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
 * groovy -cp libjar/jets3t-0.9.0.jar:libjar/httpclient-4.2-beta1.jar:libjar/httpcore-4.2-beta1.jar:libjar/commons-codec-1.4.jar:libjar/commons-logging-1.1.1.jar:libjar/java-xmlbuilder-0.4.jar upload_to_s3.groovy <file> <bucket>
 */

if( !this.args || this.args.length!=2 ){
    println 'Upload a file to S3'
    println 'usage: upload_to_s3 <file> <bucket>'
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
        largeFileObject = new S3Object(file)
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


String file = this.args[0]
String s3Out = this.args[1]
String s3Dir = ''
String bucket = s3Out

if( s3Out.indexOf('/')>=0 ){
    bucket = s3Out.substring(0, s3Out.indexOf('/'))
    s3Dir = s3Out.substring(s3Out.indexOf('/')+1)
}

if( !s3Dir.endsWith('/') )
    s3Dir += '/'
    
def finalFile = new File(file)


println "Uploading ${file}, to ${bucket+'/'+s3Dir} [${new Date()}]"
uploadLargeFile(bucket, finalFile, s3Dir+finalFile.name)

println "Finished ${new Date()}"

return 0