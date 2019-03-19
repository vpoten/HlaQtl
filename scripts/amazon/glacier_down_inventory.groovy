#!/usr/bin/env groovy

import groovy.json.JsonSlurper
import com.amazonaws.AmazonClientException
import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.policy.Policy
import com.amazonaws.auth.policy.Principal
import com.amazonaws.auth.policy.Resource
import com.amazonaws.auth.policy.Statement
import com.amazonaws.auth.policy.Statement.Effect
import com.amazonaws.auth.policy.actions.SQSActions
import com.amazonaws.services.glacier.AmazonGlacierClient
import com.amazonaws.services.glacier.model.GetJobOutputRequest
import com.amazonaws.services.glacier.model.GetJobOutputResult
import com.amazonaws.services.glacier.model.InitiateJobRequest
import com.amazonaws.services.glacier.model.InitiateJobResult
import com.amazonaws.services.glacier.model.JobParameters
import com.amazonaws.services.sns.AmazonSNSClient
import com.amazonaws.services.sns.model.CreateTopicRequest
import com.amazonaws.services.sns.model.CreateTopicResult
import com.amazonaws.services.sns.model.DeleteTopicRequest
import com.amazonaws.services.sns.model.SubscribeRequest
import com.amazonaws.services.sns.model.SubscribeResult
import com.amazonaws.services.sns.model.UnsubscribeRequest
import com.amazonaws.services.sqs.AmazonSQSClient
import com.amazonaws.services.sqs.model.CreateQueueRequest
import com.amazonaws.services.sqs.model.CreateQueueResult
import com.amazonaws.services.sqs.model.DeleteQueueRequest
import com.amazonaws.services.sqs.model.GetQueueAttributesRequest
import com.amazonaws.services.sqs.model.GetQueueAttributesResult
import com.amazonaws.services.sqs.model.Message
import com.amazonaws.services.sqs.model.ReceiveMessageRequest
import com.amazonaws.services.sqs.model.SetQueueAttributesRequest

// AWS config
String accessKey = ''
String secretKey = ''
String awsEndPoint = 'us-east-1.amazonaws.com'

if( !accessKey || !secretKey ){
    println 'Error: configure AWS credentials!'
    return 1
}

if( !this.args || this.args.length!=2 ){
    println 'Download a glacier inventory to file'
    println 'usage: glacier_down_inventory.sh <vault> <inventory_file>'
    return 1
}

String vaultName = this.args[0]
String fileName = this.args[1]


String snsTopicName = "${vaultName}_inventory_sns"
String sqsQueueName = "${vaultName}_inventory_sqs"
String sqsQueueURL
String sqsQueueARN
String snsTopicARN
String snsSubscriptionARN
long sleepTime = 120

AWSCredentials credentials = new BasicAWSCredentials(accessKey, secretKey) 

AmazonGlacierClient client = new AmazonGlacierClient(credentials)
client.setEndpoint("https://glacier." + awsEndPoint)

AmazonSQSClient sqsClient = new AmazonSQSClient(credentials)
sqsClient.setEndpoint("https://sqs." + awsEndPoint)

AmazonSNSClient snsClient = new AmazonSNSClient(credentials)
snsClient.setEndpoint("https://sns." + awsEndPoint)

//closures for AWS actions
//
def setupSQS = {
    CreateQueueRequest request = new CreateQueueRequest().withQueueName(sqsQueueName)
    CreateQueueResult result = sqsClient.createQueue(request)
    sqsQueueURL = result.getQueueUrl()

    GetQueueAttributesRequest qRequest = 
        new GetQueueAttributesRequest().withQueueUrl(sqsQueueURL).withAttributeNames("QueueArn")

    GetQueueAttributesResult qResult = sqsClient.getQueueAttributes(qRequest)
    sqsQueueARN = qResult.getAttributes().get("QueueArn")

    Policy sqsPolicy = new Policy().withStatements(
                new Statement(Effect.Allow).withPrincipals(Principal.AllUsers).
                    withActions(SQSActions.SendMessage).
                    withResources(new Resource(sqsQueueARN))
            )
    
    def queueAttributes = [:]
    queueAttributes["Policy"] = sqsPolicy.toJson()
    sqsClient.setQueueAttributes(new SetQueueAttributesRequest(sqsQueueURL, queueAttributes))

    return sqsQueueURL
}//
    

def setupSNS = {
    CreateTopicRequest request = new CreateTopicRequest().withName(snsTopicName)
    CreateTopicResult result = snsClient.createTopic(request)
    snsTopicARN = result.getTopicArn()

    SubscribeRequest request2 = 
        new SubscribeRequest().withTopicArn(snsTopicARN).withEndpoint(sqsQueueARN).withProtocol("sqs")
    SubscribeResult result2 = snsClient.subscribe(request2)

    snsSubscriptionARN = result2.getSubscriptionArn()
    
    return snsTopicARN
}//
    

def initiateJobRequest = {
        
    JobParameters jobParameters = 
        new JobParameters().withType("inventory-retrieval").withSNSTopic(snsTopicARN)

    InitiateJobRequest request = 
        new InitiateJobRequest().withVaultName(vaultName).withJobParameters(jobParameters)

    InitiateJobResult response = client.initiateJob(request)

    return response.getJobId()
}//
    

def waitForJobToComplete = { jobId, queueUrl->     
    boolean messageFound = false
    boolean jobSuccessful = false

    while (!messageFound) {
        def msgs = sqsClient.receiveMessage(
                new ReceiveMessageRequest(queueUrl).withMaxNumberOfMessages(10)
            ).messages

        msgs?.each { m->
            def slurper = new JsonSlurper()
            def jsondata = slurper.parseText( m.getBody() )

            println "Message:\n${jsondata.Message}"
            
            def json2 = slurper.parseText(jsondata.Message)
            def retrievedJobId = json2.JobId
            def statusCode = json2.StatusCode

            if ( retrievedJobId==jobId ) {
                messageFound = true
                if (statusCode == "Succeeded" ) {
                    jobSuccessful = true
                }
            }
        }
        
        Thread.sleep(sleepTime * 1000)  
    }
    return (messageFound && jobSuccessful);
}//


def downloadJobOutput = { jobId ->
        
    GetJobOutputRequest getJobOutputRequest = 
        new GetJobOutputRequest().withVaultName(vaultName).withJobId(jobId)
        
    GetJobOutputResult getJobOutputResult = client.getJobOutput(getJobOutputRequest)

    BufferedWriter out = new File(fileName).newWriter()
    BufferedReader instr = new BufferedReader(new InputStreamReader(getJobOutputResult.getBody()))         
    
    instr.eachLine{ line->
        out.write(line)
    }
    
    instr.close()
    out.close()
    
    println "Retrieved inventory to ${fileName}"
}//
    

def cleanUp = {
    snsClient.unsubscribe(new UnsubscribeRequest(snsSubscriptionARN))
    snsClient.deleteTopic(new DeleteTopicRequest(snsTopicARN))
    sqsClient.deleteQueue(new DeleteQueueRequest(sqsQueueURL))
}//


println "Start time: ${new Date()}"

try {
    def queueUrl = setupSQS()
    def snsTopic = setupSNS()
    
    println "Queue URL: ${queueUrl}"
    println "SNS topic: ${snsTopic}"

    def jobId = initiateJobRequest()
    println "Jobid = ${jobId}"

    boolean success = waitForJobToComplete(jobId, queueUrl)
    
    if (!success) {
        throw new Exception("Job did not complete successfully.")
    }

    downloadJobOutput(jobId)
}
catch (e) {
    System.err.println("Inventory retrieval failed.")
    System.err.println(e)
}
finally{
    cleanUp()
}

println "End time: ${new Date()}"

return 0
