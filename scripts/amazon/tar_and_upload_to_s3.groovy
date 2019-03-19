#!/usr/bin/env groovy

if( !this.args || this.args.length!=2 ){
    println 'Tar subfolders and upload to S3'
    println 'usage: tar_and_upload_to_s3.groovy <bucket> <directory>'
    return 1
}

String bucket = this.args[0]
String dir = this.args[1]

if( !bucket.endsWith('/') )
    bucket += '/'
    
if( !dir.endsWith('/') )
    dir += '/'
    
println "S3 bucket: ${bucket}"
println "local dir: ${dir}"

def ext = '.tar.gz'

new File(dir).eachFile{ file->
    if( file.isDirectory()  ){
        int res = "tar -czf ${file.absolutePath+ext} ${file.absolutePath}".execute().waitFor()
        
        if( res==0 ){
            def dstBucket = bucket + file.name + ext
        
            println "Copying ${file.absolutePath+ext} to ${dstBucket}"
            res = "aws put ${dstBucket} ${file.absolutePath+ext}".execute().waitFor()
            
            if( res!=0 ){
                println "Error uploading ${file.absolutePath+ext}"
            }
        }
        else{
            println "Error executing tar for ${file.absolutePath}"
        }
    }
}

return 0