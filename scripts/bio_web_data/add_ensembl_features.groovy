#!/usr/bin/env groovy

/**
 * Adds ensemble features to a TSV file using ensembl REST API
 * 
 * Calling example: add_ensembl_features.groovy ./correlation.txt.trans.filter.ext 1,4 biotype
 *
 *
 * Fields: ID, source, logic_name, external_name, feature_type, description, 
 * end, biotype, seq_region_name, strand, start, Parent
 * 
 */

import groovy.json.*

def SERVER = 'http://beta.rest.ensembl.org/'
def EXT = 'feature/id/{id}?feature=gene;feature=transcript;species=homo_sapiens;content-type=application/json'
def SUFF = '.featAdded'


if( !this.args || this.args.length!=3 ) {
    println 'Adds ensemble features to a TSV file using ensembl REST API'
    println 'Usage: add_ensembl_features.groovy <tsv file> <ids colunms> <fields>'
    return 1
}

//closures

def getEnsemblFeature = { id->
    def ext = EXT.replace('{id}',id)
    def url = (SERVER+ext).toURL()
    
    def reader = new InputStreamReader( url.openStream() )
    def slurper = new JsonSlurper()
    def result = slurper.parse(reader)
    
    return result?.find{ it['ID']==id }
}

def splitToList = { str, sep->
    return str.split(sep) as List
}

//end closures


def file = this.args[0] //input file
def columns = splitToList( this.args[1], ',').collect{it as Integer}.sort()
def fields = splitToList( this.args[2], ',')

def featMap = [:] as TreeMap
def failed = [] as TreeSet

println "Start time: ${new Date()}\n"

// 1 - retieve features
def reader = new File(file).newReader()
reader.readLine()//skip header

int count = 0

reader.splitEachLine("\t"){ toks->
    columns.each{ idx->
        def id = toks[idx]
        
        if( !(id in failed) && !featMap.containsKey(id) ){
            try{
                featMap[id] = getEnsemblFeature(id)
                count++
            } catch(e) {
                System.err.println("Error retrieving feature ${id} from ensembl REST.")
                System.err.println(e.message)
                failed << id
            }
        }
    }
}

reader.close()

println "${count} features retrieved from ensembl REST."
println "${failed.size()} features couldn't be retrieved."

// 2 - generate new file
def extHeader = columns.sum{ col-> fields.sum{"\t${it}_${col}"} }

def writer = new PrintWriter(file+SUFF)
reader = new File(file).newReader()

writer.println( reader.readLine()+extHeader )//write header

reader.eachLine{ line->
    def extLine = new StringBuilder()
    def toks = line.split("\t")
    
    columns.each{ idx->
        def feat = featMap[toks[idx]]
        fields.each{ field->
            extLine << "\t${(feat) ? feat[field] : 'NA'}" 
        }
    }
    
    writer.println( line+extLine.toString() )
}

writer.close()
reader.close()

println "End time: ${new Date()}\n"

return 0