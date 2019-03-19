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
def archiveFilter = /fastx_toolkit_(\w+)\/filtered_[12]\.fastq\.gz/ 
//example: fastx_toolkit_NA19257/filtered_1.fastq.gz

def skipSubjects = [ 
    'NA12043', 'NA12044', 'NA12144', 'NA12155', 'NA12716', 'NA12717', 'NA12750',
    'NA12813', 'NA12814', 'NA12872', 'NA12874', 'NA06994', 'NA11829', 'NA12762', 'NA12815',
    'NA06985', 'NA11830', 'NA11832', 'NA11881', 'NA11992', 'NA11993', 'NA11994', 'NA12004',
    'NA12005', 'NA12006'
] as TreeSet

def subjects = [ "NA11830",
    "NA11881",
    "NA11992",
    "NA11832",
    "NA06985",
    "NA12006",
    "NA12004",
    "NA11993",
    "NA12005",
    "NA11994",
    "NA12043",
    "NA12044",
    "NA12144",
    "NA12155",
    "NA12716",
    "NA12717",
    "NA12750",
    "NA12813",
    "NA12814",
    "NA12872",
    "NA12874",
    "NA06994",
    "NA11829",
    "NA12762",
    "NA12815",
    "NA12775",
    "NA20509",
    "NA12718",
    "NA07048",
    "NA07000",
    "NA12489",
    "NA20529",
    "NA20515",
    "NA07056",
    "NA20518",
    "NA10847",
    "NA12282",
    "NA12843",
    "NA12413",
    "NA12383",
    "NA12003",
    "NA20753",
    "NA12341",
    "NA11831",
    "NA11894",
    "NA12275",
    "NA06984",
    "NA07347",
    "NA11892",
    "NA11920",
    "NA11931",
    "NA12249",
    "NA12273",
    "NA12287",
    "NA12340",
    "NA12342",
    "NA12761",
    "NA12777",
    "NA12827",
    "NA12829",
    "NA12842",
    "NA12873",
    "NA20756",
    "NA20758",
    "NA20766",
    "NA20770",
    "NA20786",
    "NA20802",
    "NA20805",
    "NA20812",
    "HG00134",
    "HG00326",
    "HG01048",
    "HG01383",
    "NA06986",
    "NA06989",
    "NA07037",
    "NA07051",
    "NA07346",
    "NA07357",
    "NA10851",
    "NA11840",
    "NA11843",
    "NA11893",
    "NA11918",
    "NA11919",
    "NA11930",
    "NA11995",
    "NA12045",
    "NA12058",
    "NA12154",
    "NA12156",
    "NA12234",
    "NA12272",
    "NA12283",
    "NA12286",
    "NA12329",
    "NA12347",
    "NA12348",
    "NA12399",
    "NA12400",
    "NA12414",
    "NA12546",
    "NA12748",
    "NA12749",
    "NA12763",
    "NA12776",
    "NA12778",
    "NA12812",
    "NA12828",
    "NA12830",
    "NA12878",
    "NA12889",
    "NA12890",
    "NA12891",
    "NA12892",
    "NA20798",
    "NA20801",
    "NA12760"
] as TreeSet

if( !accessKey || !secretKey ){
    println 'Error: configure AWS credentials!'
    return 1
}

//call example:
//nohup ./glacier_download.sh gmsGeuvadis inventory.json fastx_toolkit_output/ > glacier_fastx.out &
if( !this.args || this.args.length!=3 ) {
    println 'Download fastx_toolkit fastq files from glacier to local'
    println 'usage: glacier_download_fastx.sh <vault> <inventory_file> <directory>'
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

def failedJobs = []


///////////////
// closures
///////////////

def download = { archiveId, downloadFilePath->
    try {
        ArchiveTransferManager atm = new ArchiveTransferManager(glacierClient, sqsClient, snsClient)
        atm.download(vaultName, archiveId, new File(downloadFilePath))
    }
    catch (Exception e) {
        println "Error downloading ${archiveId} to ${downloadFilePath}"
        return false
    }
    
    return true
}// end download closure

def testGzip = { file->
    ("gzip -t ${file}".execute().waitFor()==0 )
}// end closure


//////////////

println "Start time: ${new Date()}"

// parse JSON inventory file
def slurper = new JsonSlurper()
def jsondata = slurper.parse( new File(inventory).newReader() )

println "Vault ARN: ${jsondata.VaultARN}"
println "Inventory Date: ${jsondata.InventoryDate}"

def listFiles = jsondata.ArchiveList

println "${listFiles.size()} archives in inventory file"

if( archiveFilter ) {
    listFiles = listFiles.findAll{ 
        def mat = it.ArchiveDescription =~ archiveFilter 
        if( mat ){ return ((mat[0][1] in subjects) && !(mat[0][1] in skipSubjects)) }
        return false
    }
}

println "${listFiles.size()} archives to download after filtering"

// create out dirs if not exists
listFiles.each{ val->
    // create outdir if not exists
    def file = val.ArchiveDescription
    def dir = new File( outDir+file.substring(0, file.indexOf('/')) )
    
    if( !dir.exists() ){
        "mkdir ${dir.absolutePath}".execute().waitFor()
    }
}

// download files 
listFiles.each{ val->
    def file = val.ArchiveDescription
    
    println "Downloading ${file} : ${new Date()}"
    
    if( !download(val.ArchiveId, outDir+file) ){
        // add failed jobs to queue for repeat
        failedJobs << val
    }
    else if( file.endsWith('.gz') ){
        // check gzipped files
        if( !testGzip(outDir+file) ){
            println "Error gzip integrity for ${outDir+file}"
            failedJobs << val
        }
    }
    
}

failedJobs.each{ val->
    def file = val.ArchiveDescription
    
    println "Retry the download of ${file} : ${new Date()}"
    
    if( !download(val.ArchiveId, outDir+file) ){
        println "Error retrying the download of ${file} : ${new Date()}"
    }
    else if( file.endsWith('.gz') ){
        // check gzipped files
        if( !testGzip(outDir+file) ){
            println "Error gzip integrity for ${outDir+file}"
        }
    }
    
}


println "End time: ${new Date()}"

return 0