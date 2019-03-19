#!/usr/bin/env groovy

if( !this.args || this.args.length!=4 ){
    println 'Filters a tsv file using the selected numeric field, the field value must be less than the supplied value'
    println 'Usage: filter_by_field_numeric.groovy <src_file> <src_idx> <value> <out_file>'
    return 1
}


int srcIdx = this.args[1] as Integer
def srcFile = this.args[0]
double value = this.args[2] as Double
def outFile = this.args[3]

println "Source: ${srcFile}"
println "Field: ${srcIdx}, Value: ${value}"
println "Output: ${outFile}"

//filter source file
int total = 0
int written = 0

def reader = new File(srcFile).newReader()
def header = reader.readLine()//skip header

def writer = new File(outFile).newWriter()
writer.writeLine(header)//header

reader.eachLine{ line->
    def toks = line.split("\t")
    total++
    def current = toks[srcIdx] as Double
    
    if( current<=value /*current>=value*/){
        writer.println(line)
        written++
    }
}

reader.close()
writer.close()

println "${written} of ${total} lines written"

return 0