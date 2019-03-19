#!/usr/bin/env groovy

if( !this.args || this.args.length!=5 ){
    println 'Filters a tsv file using the selected field, the field value must be included in a list of values'
    println 'Usage: filter_by_field.groovy <src_file> <src_field_idx> <list_file> <list_field_idx> <out_file>'
    return 1
}


int srcIdx = this.args[1] as Integer
def srcFile = this.args[0]

int listIdx = this.args[3] as Integer
def listFile = this.args[2]

def outFile = this.args[4]

//get list of terms
def reader = new File(listFile).newReader()
def set = [] as TreeSet
reader.readLine()//skip header
reader.splitEachLine("\t"){ toks->
    set << toks[listIdx]
}

reader.close()

println "${set.size()} key terms read."

//filter source file
int total = 0
int written = 0

reader = new File(srcFile).newReader()
def header = reader.readLine()//skip header

def writer = new File(outFile).newWriter()
writer.writeLine(header)//header

reader.eachLine{ line->
    def toks = line.split("\t")
    total++
    
    if( toks[srcIdx] in set ){
        writer.println(line)
        written++
    }
}

reader.close()
writer.close()

println "${written} of ${total} lines written"

return 0