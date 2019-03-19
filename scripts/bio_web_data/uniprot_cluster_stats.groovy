#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

///////////////////////
// Constants
///////////////////////

def ID_100 = 'UniRef100'
def ID_90 = 'UniRef90'
def ID_50 = 'UniRef50'

def UniRefIds = [ID_100, ID_90, ID_50]
def clusters = [:]
UniRefIds.each{ clusters[it] = [:] as TreeMap }
def sizeBins = [5,10,20,30,40,50,75,100,Integer.MAX_VALUE]

def createReader = { file->
    new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(file)))))
}//

////////////////////
// Options
////////////////////

//input: uniprot id mapping file
String input = '/home/victor/work_bio/annotation/9606/HUMAN_9606_idmapping.dat.gz'
//geneontology annotation
String annotation = '/home/victor/work_bio/annotation/9606/gene_association.goa_human.gz'

// dataset generation
String clustOut = '/home/victor/Escritorio/uniprot50_clusters.txt'
double clustAnnAvg = 4.0d
int clustMax = 20
int clustMin = 10
int clUnannTol = 3


////////////////

def reader = createReader(input)

reader.splitEachLine("\t"){ toks->
    def acode = toks[0]
    def idType = toks[1]
    def id = toks[2]
    
    if( idType in UniRefIds ){
        def group = clusters[idType]
        
        def cluster = group[id]
        
        if( cluster==null ){
            cluster = []
            group[id] = cluster
        }
        
        cluster << acode
    }
}

reader.close()

clusters.each{ uniRef, group->
    println "\n===== ${uniRef}: ${group.size()} clusters"
    
    def counts = [:]
    sizeBins.each{ counts[it]=0 }
    
    group.each{ id, cluster->
        def bin = sizeBins.find{ cluster.size()<=it }
        counts[bin] = counts[bin]+1
    }
    
    sizeBins.each{ println "${counts[it]} clusters <= ${it} members" }
}

// parse annotation goa file
def indexGOA = [:] as TreeMap
reader = createReader(annotation)

reader.eachLine{ line->
    if( !line.startsWith('!') ){
        def toks = line.split("\t",-1)
        def acode = toks[1]
        def go = toks[4]

        def list = indexGOA[acode]

        if( list==null ){
            list = []
            indexGOA[acode] = list
        }

        list << go
    }
}

reader.close()

// annotation stats for UniRef50 clusters
println "\n===== UniRef50 geneontology annotation:"

def counts = [:]
def annCounts = [:]
sizeBins.each{ counts[it]=0; annCounts[it]=0.0d }

def selectedCl = [] as TreeSet

clusters[ID_50].each{ id, cluster->
    int unannotated = cluster.sum{ac-> (indexGOA[ac]) ? 0 : 1}
    
    double num = cluster.sum{ ac-> 
        def annot = indexGOA[ac] 
        (annot) ? annot.size() : 0.0
    }
    
    int size = cluster.size()
    num /= (double)size
    
    if( unannotated<=clUnannTol && num>=clustAnnAvg && size>=clustMin && size<=clustMax ){
       selectedCl << id 
    }

    def bin = sizeBins.find{ cluster.size()<=it }
    counts[bin] = counts[bin]+1
    annCounts[bin] = annCounts[bin] + num
}

sizeBins.each{
    println "${annCounts[it]/(double)counts[it]} average annnotations in clusters <= ${it} members (${counts[it]})"
}

println "${selectedCl.size()} clusters meet selection criteria: [${clustMin},${clustMax}], annAvg>=${clustAnnAvg}"


if( clustOut ){
    println "\nWriting clusters to: ${clustOut}"
    def writer = new PrintWriter(clustOut)
    
    clusters[ID_50].each{ id, cluster->
        if( id in selectedCl){
            cluster.each{ ac->
                if( indexGOA[ac] ){
                    writer.println("${id}\t${ac}\t${indexGOA[ac].sum{it+','}}")
                }
            }
        }
    }
    
    writer.close()
}

return 0