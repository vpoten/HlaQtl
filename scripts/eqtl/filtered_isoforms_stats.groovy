#!/usr/bin/env groovy


int F_ISO_ID = 0
int F_ISO_GENE = 5
int F_ISO_TYPE = 7

def BIOTYPE_TRANS = ['IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IG_LV_gene',
        'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay',
        'polymorphic', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene','TR_J_gene', 
        'TR_V_gene', 'miRNA', 'antisense', 'antisense_RNA'] as TreeSet


if( !this.args || this.args.length!=2 ){
    println 'Calculates stats and filters by biotype a filtered+ensembl_feat isoform file.'
    println 'Usage: filtered_isoforms_stats.groovy <file> <out_file>'
    return 1
}

def input = this.args[0]
def output = this.args[1]


println "Input file: ${input}"
println "Output file: ${output}"

def genes = [:] as TreeMap
int count = 0

def reader = new File(input).newReader()
def header = reader.readLine()//skip header

reader.eachLine{ line->
    def toks = line.split("\\s")
    def gene = toks[F_ISO_GENE]
    
    if( gene!='NA' && (toks[F_ISO_TYPE] in BIOTYPE_TRANS) ) {
        if( !genes.containsKey(gene) ){
            genes[gene] = []
        }
        
        genes[gene] << line
        count++
    }
}

reader.close()

// print out results
println "${genes.size()} genes after final filtering."
println "${count} transcripts after final filtering."

def writer = new File(output).newWriter()
writer.writeLine(header)//header

genes.each{ gene, list->
    list.each{ writer.writeLine(it) }
}

writer.close()

return 0