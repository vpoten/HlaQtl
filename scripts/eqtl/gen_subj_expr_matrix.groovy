#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

String MERG_TRACK_FILE = 'iso_merged.tracking'
String TFAM_FILE = '1000g_snps.tfam'

//merged_tracking file fields: subj iso iso_merg locus fpkm
int F_SUBJ = 0
int F_ISO = 2
int F_FPKM = 4


//closures

def getFileIfCompress = { name->
    def file = new File(name)

    if( !file.exists() )
        file = new File(name+'.gz')
    else
        return file

    if( !file.exists() )
        return null

    return file
}

def createReader = { file->
    (file.name.endsWith('.gz')) ? 
        new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))) : file.newReader()
}

def subjectList = { dir->
    def dirs = []
    new File(dir).eachDir{ dirs << it }
    
    def tfamDir = dirs.find{ (new File(it.absolutePath+'/'+TFAM_FILE)).exists() }

    def list = []

    // read subject id from .tfam file
    def tfamReader = new File(tfamDir.absolutePath+'/'+TFAM_FILE).newReader()
    tfamReader.splitEachLine("\\s"){ toks-> list << toks[0] }
    tfamReader.close()

    return list
}
    
//end closures


if( !this.args || this.args.length!=3 ) {
    println "Generates subjects expression matrix file: GxN, G=subjects, N=isoforms expr."
    println 'Usage: get_subj_expr_matrix.groovy <eqtl result dir> <isoforms file> <out file>'
    return 1
}

def inputdir = this.args[0]
def isoFile = this.args[1]
def outfile = this.args[2]


println "Input dir: ${inputdir}"
println "Outfile: ${outfile}"
println "Filtered isoforms: ${isoFile}"

println "Start time: ${new Date()}\n"

//parse selected isoforms filter file
def selectedIso = [] as TreeSet

def reader = new File(isoFile).newReader()
reader.readLine()//skip header
reader.splitEachLine("\t"){ toks-> selectedIso << toks[0] }
reader.close()

println "${selectedIso.size()} isoforms read"

//load subjects
def subjects = subjectList(inputdir)
println "${subjects.size()} subjects read"

def isoExprData = [:] as TreeMap //key=subjectId, value={map with key=iso,value=expresion}
subjects.each{ isoExprData[it] = [:] as TreeMap }

// traverse result subfolders
new File(inputdir).eachDir{ dir->
    File mergTrackFile = getFileIfCompress(dir.absolutePath+'/'+MERG_TRACK_FILE)

    if( mergTrackFile ){
        reader = createReader(mergTrackFile)

        reader.splitEachLine("\\s"){ toks->
            //file fields: subj iso iso_merg locus fpkm
            if( (toks[F_ISO] in selectedIso) && isoExprData.containsKey(toks[F_SUBJ]) ){
                isoExprData[toks[F_SUBJ]][toks[F_ISO]] = toks[F_FPKM] as Double
            }
        }

        reader.close()
    }
}

// write expression matrix
def writer = new PrintWriter(outfile)
int end = selectedIso.size()-1

//print header
writer.print('sample,')
selectedIso.eachWithIndex{ iso, i->
    writer.print( iso+((i==end) ? '\n' : ',') )
}

isoExprData.each{ subj, isoforms->
    writer.print(subj+',')//print row name (subject id)
    
    selectedIso.eachWithIndex{ iso, i->
        writer.print( (isoforms[iso] ?: 0.0)+((i==end) ? '\n' : ',') )
    }
}
writer.close()

println "\nEnd time: ${new Date()}\n"

return 0