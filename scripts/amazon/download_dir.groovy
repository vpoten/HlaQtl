#!/usr/bin/env groovy

if( !this.args || this.args.length!=2 ){
    println 'Download a directory from S3 to local'
    println 'usage: download_dir <S3_input_dir> <local_output_dir>'
    return 1
}

String dir = this.args[0]
String outdir = this.args[1]

if( !dir.endsWith('/') )
    dir += '/'

if( !outdir.endsWith('/') )
    outdir += '/'
    
println "S3 dir: ${dir}"
println "Output dir: ${outdir}"
    
def proc = "aws ls -1 ${dir}".execute()

if( proc.waitFor()!=0 ){
    System.err.println("Cannot list ${dir}")
    return 1
}
    
new StringReader(proc.text).eachLine{ line->
    def name = line.substring(line.lastIndexOf('/')+1)
    println "downloading ${dir+name}"
    if( "aws get ${dir+name} ${outdir}".execute().waitFor()!=0 ){
        System.err.println("Cannot download ${dir+name}")
    }
}

return 0
