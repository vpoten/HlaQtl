#!/usr/bin/env groovy

String CORREL_FILE = 'correlation.txt.filter'
double pValThr = 1e-5

//correlation.txt file fields
int F_SNP = 0
int F_ISO = 1
int F_CORR = 2
int F_PVAL = 5


if( !this.args || this.args.length!=2 ) {
    println "Get best eqtls from a eqtl results dir."
    println "Parses correlation.txt.filter files to get the best snps per isoform."
    println 'Usage: get_best_eqtl.groovy <eqtl result dir> <out file pref>'
    return 1
}

def inputdir = this.args[0]
def outfile = this.args[1]

println "Input dir: ${inputdir}"
println "Outfile: ${outfile}"
println "p-value threshold: ${pValThr}"

println "Start time: ${new Date()}\n"

def snpsSet = [] as TreeSet

def parseCorrFile = { file, chr->
    def reader = file.newReader()
    reader.readLine()//skip header
    
    def isoCorrMap = [:]

    reader.splitEachLine("\\s"){ toks->
        double pval = toks[F_PVAL] as Double
        
        if( pval<pValThr ){
            //store snp-iso correlation value
            if( isoCorrMap[toks[F_ISO]]==null )
                isoCorrMap[toks[F_ISO]] = []
                
            isoCorrMap[toks[F_ISO]] << [ (F_SNP):toks[F_SNP],
                (F_CORR):(toks[F_CORR] as Double), (F_PVAL):pval ]
        }
    }

    reader.close()
    
    //get the best snp of each isoform
    isoCorrMap.each{ iso, list->
        def bestSnp = list.max{ Math.abs(it[F_CORR]) }
        snpsSet << "${chr}\t${bestSnp[F_SNP]}"
        println "${iso}\t${bestSnp[F_SNP]}\t${bestSnp[F_CORR]}\t${bestSnp[F_PVAL]}\t${chr}"
    }
    
}// end of parseCorrFile

// traverse result subfolders
new File(inputdir).eachDir{ dir->
    def file = new File(dir.absolutePath+'/'+CORREL_FILE)
    
    if( file.exists() ){
        def chr = dir.name.substring(0, dir.name.indexOf('_'))
        parseCorrFile(file, chr)
    }
}

// write filtered snps
writer = new PrintWriter(outfile)
snpsSet.each{ writer.println(it) }
writer.close()

println "\nEnd time: ${new Date()}\n"

return 0