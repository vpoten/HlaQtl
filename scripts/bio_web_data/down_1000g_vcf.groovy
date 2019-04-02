#!/usr/bin/env groovy

String _1000G_FTP = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
String S3_VCF_FILE = 'ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
String S3_VCF_TBI = 'ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi'
String WGET_COMMAND = "wget -nv --no-proxy"

// def chrs = (1..22).collect{ 'chr'+it }
def chrs = [1, 3, 5, 7].collect{ 'chr'+it }
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