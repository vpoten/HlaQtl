#!/usr/bin/env groovy

import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.services.glacier.AmazonGlacierClient
import com.amazonaws.services.glacier.transfer.ArchiveTransferManager
import com.amazonaws.services.glacier.transfer.UploadResult

// AWS config
String accessKey = ''
String secretKey = ''
String awsEndPoint = 'https://glacier.us-east-1.amazonaws.com/'

if( !accessKey || !secretKey ){
    println 'Error: configure AWS credentials!'
    return 1
}

if( !this.args || this.args.length!=2 ){
    println 'Upload a directory from local to glacier'
    println 'usage: glacier_upload.sh <vault> <directory>'
    return 1
}



String vault = this.args[0]
String dir = this.args[1]

if( !dir.endsWith('/') ){ dir += '/' }
    
println "Glacier vault: ${vault}"
println "Local dir: ${dir}"

// AWS objects
AmazonGlacierClient client
AWSCredentials credentials = new BasicAWSCredentials(accessKey, secretKey) 

client = new AmazonGlacierClient(credentials)
client.setEndpoint(awsEndPoint)

//upload file to glacier
def uploadFile = { file->
    try {
        ArchiveTransferManager atm = new ArchiveTransferManager(client, credentials)
        String desc = file.absolutePath.substring(dir.length())
        println "Uploading ${file.absolutePath} to vault ${vault}"
        UploadResult result = atm.upload(vault, desc, file)
        println "Archive ID for ${desc}: ${result.getArchiveId()}"

    } catch (e) {
        println "Error uploading ${file.absolutePath}"
        System.err.println(e);
    }
}


println "Start time: ${new Date()}"

new File(dir).eachFileRecurse{ file->
    if( file.file )
        uploadFile(file)
}

println "End time: ${new Date()}"

return 0