#!/usr/bin/env groovy

import groovy.json.JsonSlurper
import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.services.glacier.AmazonGlacierClient
import com.amazonaws.services.glacier.transfer.ArchiveTransferManager
import com.amazonaws.services.sns.AmazonSNSClient
import com.amazonaws.services.sqs.AmazonSQSClient

// AWS config
String accessKey = ''
String secretKey = ''
String awsEndPoint = 'us-east-1.amazonaws.com'
String archiveFilter = '.bam'

if( !accessKey || !secretKey ){
    println 'Error: configure AWS credentials!'
    return 1
}

if( !this.args || this.args.length!=3 ){
    println 'Upload an inventory from glacier to local'
    println 'usage: glacier_download.sh <vault> <inventory_file> <directory>'
    return 1
}

String vaultName = this.args[0]
String inventory = this.args[1]
String outDir = this.args[2]
   
if( !outDir.endsWith('/') ){ outDir += '/' } 
        
println "Glacier vault: ${vaultName}"
println "Local dir: ${outDir}"
println "Inventory file: ${inventory}"

AWSCredentials credentials = new BasicAWSCredentials(accessKey, secretKey) 

AmazonGlacierClient glacierClient = new AmazonGlacierClient(credentials)

AmazonSQSClient sqsClient = new AmazonSQSClient(credentials)
AmazonSNSClient snsClient = new AmazonSNSClient(credentials)
glacierClient.setEndpoint("glacier."+awsEndPoint)
sqsClient.setEndpoint("sqs."+awsEndPoint)
snsClient.setEndpoint("sns."+awsEndPoint)


def download = { archiveId, downloadFilePath->
    try {
        ArchiveTransferManager atm = new ArchiveTransferManager(glacierClient, sqsClient, snsClient)
        atm.download(vaultName, archiveId, new File(downloadFilePath))
    }
    catch (Exception e) {
        println "Error downloading ${archiveId} to ${downloadFilePath}"
        System.err.println(e)
    }
}// end download closure


println "Start time: ${new Date()}"

// parse JSON inventory file
def slurper = new JsonSlurper()
def jsondata = slurper.parse( new File(inventory).newReader() )

println "Vault ARN: ${jsondata.VaultARN}"
println "Inventory Date: ${jsondata.InventoryDate}"

def listFiles = jsondata.ArchiveList

println "${listFiles.size()} archives in inventory file"

if( archiveFilter ){
    listFiles = listFiles.findAll{ it.ArchiveDescription.contains(archiveFilter) }
}

println "${listFiles.size()} archives to download after filtering"

listFiles.each{ val->
    // create outdir if not exists
    def file = val.ArchiveDescription
    def dir = new File( outDir+file.substring(0, file.indexOf('/')) )
    
    if( !dir.exists() ){
        "mkdir ${dir.absolutePath}".execute().waitFor()
    }
    
    println "Downloading ${file} : ${new Date()}"
    download(val.ArchiveId, outDir+file)
}

println "End time: ${new Date()}"

return 0