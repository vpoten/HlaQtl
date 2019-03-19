#!/usr/bin/env groovy

/**
 * Adds ensemble features to a TSV file using ensembl REST API
 * 
 * Calling example: repeat_ensembl_features.groovy file_ensembl.out Parent,external_name,biotype,seq_region_name,start,end
 *
 *
 * Fields: ID, source, logic_name, external_name, feature_type, description, 
 * end, biotype, seq_region_name, strand, start, Parent
 * 
 */

import groovy.json.*

def SERVER = 'http://beta.rest.ensembl.org/'
def EXT = 'feature/id/{id}?feature=gene;feature=transcript;species=homo_sapiens;content-type=application/json'


if( !this.args || this.args.length!=2 ) {
    println 'Repeat ensembl features retrieval for timed ou'
    println 'Usage: repeat_ensembl_features.groovy <error file> <fields>'
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
def fields = splitToList( this.args[1], ',')


// 1 - retieve features
def reader = new File(file).newReader()
def prevLine = reader.readLine()

reader.eachLine{ line->
    if( line.contains('timed out') ){
        def id = prevLine.substring( prevLine.indexOf('ENST'), prevLine.indexOf('from')-1 )

        try{
            def map = getEnsemblFeature(id)
            def extLine = new StringBuilder()
            extLine << "${id}"
            fields.each{ field->
                extLine << "\t${(map) ? map[field] : 'NA'}" 
            }

            println extLine.toString()
        } catch(e) {
            System.err.println("Error retrieving feature ${id} from ensembl REST: ${e.message}")
        }
    }
    
    prevLine = line
}

reader.close()


return 0