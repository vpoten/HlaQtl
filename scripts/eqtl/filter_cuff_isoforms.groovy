#!/usr/bin/env groovy

import groovy.json.*
import java.util.zip.GZIPInputStream

// Ensembl REST feature fields: ID, source, logic_name, external_name, feature_type, description, 
// end, biotype, seq_region_name, strand, start, Parent
def SERVER = 'http://beta.rest.ensembl.org/'
def EXT = 'feature/id/{id}?feature=gene;feature=transcript;species=homo_sapiens;content-type=application/json'


int F_ISO_ID = 0
int F_ISO_FPKM = 9
int F_ISO_COVER = 8
int F_STATUS = 12

def BIOTYPE_TRANS = ['IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IG_LV_gene',
        'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay',
        'polymorphic', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene','TR_J_gene', 
        'TR_V_gene', 'miRNA', 'antisense', 'antisense_RNA'] as TreeSet

float fpkmThr = 3.0
double freqThr = 0.1

def regexCuffOut = /(\w+)\.fpkm_tracking[\.gz]*/

//closures

def getEnsemblFeature = { id->
    def ext = EXT.replace('{id}',id)
    def url = (SERVER+ext).toURL()
    
    def reader = new InputStreamReader( url.openStream() )
    def slurper = new JsonSlurper()
    def result = slurper.parse(reader)
    
    return result?.find{ it['ID']==id }
}

//end closures

if( !this.args || this.args.length!=2 ){
    println 'Filter cufflinks isoforms getting real transcript isoforms.'
    println 'Usage: filter_cuff_isoforms.groovy <data dir> <out_file>'
    return 1
}

def input = this.args[0]
def output = this.args[1]

if( !input.endsWith('/') ){ input += '/' }

println "Input dir: ${input}"
println "Output file: ${output}"
println "Frequency thr.: ${freqThr}"
println "FPKM thr.: ${fpkmThr}"

println "Start time: ${new Date()}"


// get FPKM of each isoform
def files = new File(input).list({d, f-> f ==~ regexCuffOut } as FilenameFilter).toList()

def fpkms = [:] as TreeMap

println "${files.size()} FPKM cufflinks files to process."

files.each{ name->
    def file = new File(input+name)
    def reader = (name.endsWith('.gz')) ? 
        new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
    
    reader.readLine()//skip header
    
    reader.splitEachLine("\\s"){ toks->
        if( toks[F_STATUS]=='OK' ){
            def id = toks[F_ISO_ID]
            
            if( !fpkms.containsKey(id) ){
                fpkms[id] = []
            }
            
            fpkms[id] << (toks[F_ISO_FPKM] as Float)
        }
    }
    
    reader.close()
}

def selected = [] as Set

// get isoforms that satisfy the frequency/FPKM threshold
fpkms.each{ id, list->
    double count = list.sum{ (it>=fpkmThr) ? 1.0 : 0.0 }
    double rate = count / ((double) list.size())
    
    if( rate >= freqThr ){
        selected << id
    }
}

println "${selected.size()} isoforms selected after frequency filtering."


// print out results

def writer = new File(output).newWriter()
writer.writeLine("isoform\texpr_avg\texpr_max\texpr_min\t#subjs")//header

selected.each{ id->
    def list = fpkms[id]
    count = list.sum{ (it>=fpkmThr) ? 1 : 0 }
    double avg = (list.sum()/(double)list.size())
    
    writer.writeLine("${id}\t${avg}\t${list.max()}\t${list.min()}\t${count}")
}

writer.close()

println "End time: ${new Date()}"

return 0