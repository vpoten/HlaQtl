#!/usr/bin/env groovy

if( !this.args || this.args.length!=2 ){
    println 'Transform bedGraph'
    println 'Usage: transform_bedgraph.groovy <file in bedGraph> <file out bedGraph>'
    return 1
}

def fileIn = this.args[0]
def fileOut = this.args[1]
int width = 1

println "File input: ${fileIn}"
println "File output: ${fileOut}"

def reader = new File(fileIn).newReader()
def writer = new File(fileOut).newPrintWriter()

writer.println( reader.readLine() )//write header

def currChr = null
def regions = []

def flushData = {
    int start = regions[0][0]
    int end = regions[0][0]+width
    double score = regions[0][1]
    
    for(int i=1; i<regions.size(); i++){
        
        if( end >= regions[i][0] ){
            end = regions[i][0]+width
            score = Math.max(score, regions[i][1])
            
        }
        else{
            writer.println("${currChr}\t${start}\t${end}\t${score}")
            start = regions[i][0]
            end = regions[i][0]+width
            score = regions[i][1]
        }
    }
    
    if( regions.size()>1 )
        writer.println("${currChr}\t${start}\t${end}\t${score}")
    
    regions.clear()
}//end flushData closure

reader.splitEachLine("\\s"){ toks->
    
    if( currChr!=toks[0] ){
        if( currChr!=null ){
            flushData()
        }
        currChr = toks[0]
    }
    regions << [toks[1] as Integer, toks[3] as Double]
}

flushData()
    
reader.close()
writer.close()

return 0
