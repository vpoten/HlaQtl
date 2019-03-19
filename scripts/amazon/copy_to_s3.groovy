#!/usr/bin/env groovy

if( !this.args || this.args.length!=2 ){
    println 'Copy a directory from local to S3 recursively'
    println 'usage: copy_to_s3.sh <bucket> <directory>'
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


new File(dir).eachFileRecurse{ file->
    if( file.file ){
        def dstBucket = bucket + file.absolutePath.substring(dir.length())
        println "Copying ${file.absolutePath} to ${dstBucket}"
        "aws put ${dstBucket} ${file.absolutePath}".execute().waitFor()
    }
}

return 0