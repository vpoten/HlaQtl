#!/usr/bin/env groovy

String _1000G_FTP = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/'
String S3_VCF_FILE = 'ALL.{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz'
String S3_VCF_TBI = 'ALL.{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi'
String WGET_COMMAND = "wget -nv --no-proxy"

def chrs = (1..22).collect{ 'chr'+it }
chrs << 'chrX'

String outdir = this.args[0]

if( !outdir.endsWith('/') )
    outdir += '/'
    
// download 100g file
def downUrl = { file, chr->
    file = file.replace('{chr}',chr)
    println "downloading ${file}"
    def url = _1000G_FTP+file
    
    if( "${WGET_COMMAND} -O ${outdir+file} ${url}".execute().waitFor()!=0 )
        println "Error downloading ${file}"
}

chrs.each{ chr->
    [S3_VCF_FILE,S3_VCF_TBI].each{ file->
        downUrl(file,chr)
    }
}

return 0