/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import org.awsjoblauncher.storage.*


/**
 * helper class
 * 
 * @author victor
 */
class Genotype {
    String subject //subject id
    def snps = [:] as TreeMap //map with key=SNP_id and value=alleles
}

/**
 *
 * @author victor
 */
class SNPManager {
    //constants
    public static final String filtOpts =  '--hwe-all --maf 0.05 --hwe 0.05' //plink filtering options
    
    static final String HAPMAP_URL = 'http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/latest/forward/non-redundant/'
    static final String hapmap_file = "genotypes_{chr}_{group}_r27_nr.b36_fwd.txt.gz"
    static final def HAPMAP_GROUPS = ['CEU','ASW','CHB',/*'CHD','GIH','JPT','LWK',*/'MEX',/*'MKK','TSI',*/'YRI']
    
    static final String HAPMAP_PHASED_URL = 'ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2009-02_phaseIII/HapMap3_r2/{group}/'
    static final String trios_phased_file = 'hapmap3_r2_b36_fwd.consensus.qc.poly.{chr}_{lgroup}.phased.gz'
    static final String duos_phased_file = 'hapmap3_r2_b36_fwd.consensus.qc.poly.{chr}_{lgroup}.D.phased.gz'
    static final String unrel_phased_file = 'hapmap3_r2_b36_fwd.consensus.qc.poly.{chr}_{lgroup}.unr.phased.gz'
    
    static final String HAPMAP_PED = 'hapmap_snps'
    static final String HAPMAP_HAP = 'hapmap_haplo'
    static final String _1000G_PED = '1000g_snps'
    static final String _1000G_HAP = '1000g_haplo'
    static final String MERGED_PED = 'merged'
    static final String MERGED_HAP = 'merged_hap'
    static final String MERGED_FIN = 'merged_final'
    static final String IMPUTED_PRE = 'imputed'
    static final String IMPUTED_FIN = 'imputed_final'
    static final String SLICE_VCF = 'slice.vcf'
    static final String BEAGLE_JAR = 'beagle.jar'
    static final double R2_CUTOFF = 0.5
    
    //plink executable; overwrite if necessary
    static String PLINK = 'plink'
    static String VCFTOOLS = 'vcftools'
    static String TABIX = 'tabix'
    
    /**
    HapMap File format.
    In general, genotype files consists of tab-delimited text 
    files with at least the following columns:

    Col1: refSNP rs# identifier at the time of release (NB might merge 
          with another rs# in the future)
    Col2: SNP alleles according to dbSNP
    Col3: chromosome that SNP maps to 
    Col4: chromosome position of SNP, in basepairs on reference sequence
    Col5: strand of reference sequence that SNP maps to
    Col6: version of reference sequence assembly
    Col7: HapMap genotype center that produced the genotypes
    Col8: LSID for HapMap protocol used for genotyping
    Col9: LSID for HapMap assay used for genotyping
    Col10: LSID for panel of individuals genotyped
    Col11: QC-code, currently 'QC+' for all entries (for future use)
    Col12 and on: observed genotypes of samples, one per column, sample 
    identifiers in column headers (Coriell catalog numbers, example: 
    NA10847). Duplicate samples have .dup suffix.
     */
    
    /**
     The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype (-9 missing;0 missing;1 unaffected;2 affected)
     
     Map file:
     chromosome (1-22, X, Y or 0 if unplaced)
     rs# or snp identifier
     Genetic distance (morgans)
     Base-pair position (bp units)
     */
    
    static final String _1000G_FTP = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/'
    static final String S3_VCF_FILE = 'ALL.{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz'
    static final String S3_VCF_TBI = 'ALL.{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi'
    static final String S3_VCF_URL = 's3://1000genomes/release/20110521/'
    static def regexVcf = /ALL\.(\w+)\.phase1_release_v3\.20101123\.snps_indels_svs\.genotypes\.vcf\.gz/
    
    def hmSnps = [:] as TreeMap //map of selected SNPs {snp_id, SNPData} (HapMap)
    def hmSortSnps = [] //snps sorted by position (HapMap)
    def genotypes = [:] as TreeMap //map of subjects genotypes {subject_id,Genotype}
    
    def hmHaplo = [:] as TreeMap //map of selected SNPs {snp_id, SNPData} (HapMap phased)
    def hmSortHaplo = [] //snps sorted by position (HapMap phased)
    def haplotypes = [:] as TreeMap //map of subjects genotypes {subject_id,haplotype}
    
    def imputedSnps = [:] as TreeMap //map of merged SNPs {snp_id, SNPData}
    def imputedSort = [] // merged snps sorted by position
    def discardedSnps = [] as TreeSet //set of not reliable imputed snps
    String workDir
    String chr //chromosome
    
    
    /**
     * @param subjects : a list of subject ids
     * @return a SNPManager object
     */
    static loadSNPData(subjects, outdir, vcfFile, List groups, String locusStr, boolean _1000gOnly, snpsFile = null) {
        def locus = ((locusStr!=null) ? (locusStr =~ Utils.locusRegex) : null)
        def chr = (locusStr!=null) ? locus[0][1] : null
        def start = (locusStr!=null) ? (locus[0][2] as Integer) : null
        def end = (locusStr!=null) ? (locus[0][3] as Integer) : null
        
        def snpMng = new SNPManager(workDir:outdir, chr:chr)
        
        String hmgroup = null
        
        if( groups.size()==1 ){
            if( groups[0] in HAPMAP_GROUPS )
                hmgroup = groups[0]
        }
        
        
        // download 1000 genomes data
        def tbiFile = null
        def tmpDir = outdir+'tmp/'
        Utils.createDir(tmpDir)
        
        if( !vcfFile ){
            // get vcf from 1000genomes ftp using tabix
            vcfFile = _1000G_FTP+S3_VCF_FILE
            tbiFile = _1000G_FTP+S3_VCF_TBI
            vcfFile = vcfFile.replace('{chr}',chr)
            tbiFile = tbiFile.replace('{chr}',chr)
        }
        
        if( snpsFile && locusStr==null ){
            //trans-eqtl calculations (vcfFile is the dir where .vcf are stored)
            if( !vcfFile.endsWith('/') )
                vcfFile += '/'
                
            snpMng.vcfDirToTpedAndBgl(vcfFile, outdir+_1000G_PED, subjects, snpsFile)
        }
        else if( !_1000gOnly && hmgroup ) {
            //download and parse hapMap genotype data
            def file = snpMng.downloadSNPs(hmgroup)

            if( !file ){
                System.err.println("Error downloading HapMap genotypes")
                return null
            }

            snpMng.parseHapMap(file, start, end)

            //download and parse hapMap phasing data
            if( !snpMng.downloadPhasingSNPs(hmgroup) ){
                System.err.println("Error downloading HapMap phasing")
                return null
            }

            snpMng.parseHapMapPhasing(start, end, hmgroup)

            //remove downloaded files
            snpMng.clearHapMapFiles(hmgroup)
            
            //count subjects present in hapMap
            int num = subjects.sum{ (it in snpMng.genotypes.keySet()) ? 1 : 0 }
            def missing = subjects.findAll{ !(it in snpMng.genotypes.keySet()) }
            println "${num} subjects of ${subjects.size()} present in HapMap set of ${snpMng.genotypes.keySet().size()} individuals"
            println "missing: ${missing}"
            
            // convert to ped
            snpMng.hapMapToPed( outdir+HAPMAP_PED, snpMng.genotypes, snpMng.hmSortSnps)
            
            // convert to ped
            vcfToPed( vcfFile, tbiFile, chr, start, end, outdir+_1000G_PED, false,
                snpMng.genotypes.keySet()+subjects, tmpDir )

            // generate 1000genomes phased haplotype ped/map files
            vcfToPed( vcfFile, tbiFile, chr, start, end, outdir+_1000G_HAP, true, 
                snpMng.genotypes.keySet()+subjects, tmpDir )
        
            //merge and filter ped files
            if( !mergePedFiles(outdir+HAPMAP_PED, outdir+_1000G_PED, outdir+MERGED_PED) )
                return null

            snpMng.generateBgl(outdir+MERGED_PED)
        
            // generate HapMap phased ped/map files
            snpMng.hapMapToPed( outdir+HAPMAP_HAP, snpMng.haplotypes, snpMng.hmSortHaplo)

            // merge 1000genomes and HapMap phased
            if( !mergePedFiles(outdir+HAPMAP_HAP, outdir+_1000G_HAP, outdir+MERGED_HAP) )
                return null

            snpMng.generatePhasedBgl(outdir+MERGED_HAP)
        
            // impute missing genotypes using phased data
            def beagleJar = installBeagle(tmpDir)
            runBeagle( beagleJar, ["unphased=${outdir+MERGED_PED+'.bgl'}",
                "phased=${outdir+MERGED_HAP+'.bgl'}","markers=${outdir+MERGED_HAP+'.mark'}",
                'missing=0',"out=${outdir+IMPUTED_PRE}"] )

            // parse and filter first result
            snpMng.parseBgl("${outdir}${IMPUTED_PRE}.${MERGED_PED}.bgl.phased.gz", R2_CUTOFF)

            // add not phasing SNPs to imputation result
            snpMng.generateMergedBgl(outdir+MERGED_FIN, outdir+MERGED_PED)

            // impute missing genotypes using unphased data
            runBeagle( beagleJar, ["unphased=${outdir+MERGED_FIN+'.bgl'}",
                'missing=0',"out=${outdir+IMPUTED_FIN}"] )

            // parse and filter final result
            snpMng.parseBgl("${outdir}${IMPUTED_FIN}.${MERGED_FIN}.bgl.phased.gz", R2_CUTOFF)
        }
        else {
            // use 1000 genomes only
            
            // convert to tped
            vcfToTped( vcfFile, tbiFile, chr, start, end, outdir+_1000G_PED, false,
                subjects, tmpDir )
            
            if( !filterTpedFiles(outdir+_1000G_PED) )
                return null
               
            // generate bgl file and populates imputedSnps map
            snpMng.generateBgl(outdir+_1000G_PED) 
        }
        
        // generate file with imputed filtered results
        snpMng.writeGenotypes("${outdir}${IMPUTED_FIN}.filter.txt")
            
        //delete tmp dir
        "rm -rf ${tmpDir}".execute().waitFor()
        
        return snpMng
    }
    
    /**
     * 
     */ 
    static loadHapMapSNPData(subjects, outdir, List groups, String locusStr) {
        def locus = ((locusStr!=null) ? (locusStr =~ Utils.locusRegex) : null)
        def chr = (locusStr!=null) ? locus[0][1] : null
        def start = (locusStr!=null) ? (locus[0][2] as Integer) : null
        def end = (locusStr!=null) ? (locus[0][3] as Integer) : null
        
        def snpMng = new SNPManager(workDir:outdir, chr:chr)
        
        String hmgroup = null
        
        if( groups.size()==1 ){
            if( groups[0] in HAPMAP_GROUPS )
                hmgroup = groups[0]
        }
        
        def tmpDir = outdir+'tmp/'
        Utils.createDir(tmpDir)
        
        //download and parse hapMap genotype data
        def file = snpMng.downloadSNPs(hmgroup)

        if( !file ){
            System.err.println("Error downloading HapMap genotypes")
            return null
        }

        snpMng.parseHapMap(file, start, end)
        
        //remove downloaded files
        snpMng.clearHapMapFiles(hmgroup)

        //count subjects present in hapMap
        int num = subjects.sum{ (it in snpMng.genotypes.keySet()) ? 1 : 0 }
        def missing = subjects.findAll{ !(it in snpMng.genotypes.keySet()) }
        println "${num} subjects of ${subjects.size()} present in HapMap set of ${snpMng.genotypes.keySet().size()} individuals"
        println "missing: ${missing}"

        // convert to ped
        snpMng.hapMapToPed( outdir+HAPMAP_PED, snpMng.genotypes, snpMng.hmSortSnps)

        //filter and covert to tped
        def test = "--recode ${filtOpts} --transpose "
        def comm = plinkTest(test, outdir+HAPMAP_PED+'.ped', outdir+HAPMAP_PED+'.map', outdir+HAPMAP_PED)
        
        if( comm.execute().waitFor()!=0 ){
            System.err.println("Error recoding plink files [${outdir+HAPMAP_PED}]")
            return null
        }

        snpMng.generateBgl(outdir+HAPMAP_PED)
        
        // generate file with imputed filtered results
        snpMng.writeGenotypes("${outdir}${IMPUTED_FIN}.filter.txt")
            
        //delete tmp dir
        "rm -rf ${tmpDir}".execute().waitFor()
        
        return snpMng
    }
    
    /**
     *
     */
    protected static def mergePedFiles(file1Base, file2Base, outBase){
        def test = "--merge ${file2Base+'.ped'} ${file2Base+'.map'} --merge-mode 2 --recode ${filtOpts}"
        def comm = plinkTest(test, file1Base+'.ped', file1Base+'.map', outBase)
        
        if( comm.execute().waitFor()!=0 ){
            System.err.println("Error merging/recoding plink files")
            return null
        }
        
        test = "--recode --transpose"
        comm = plinkTest(test, outBase+'.ped', outBase+'.map', outBase)
        
        if( comm.execute().waitFor()!=0 ){
            System.err.println("Error merging/recoding plink files")
            return null
        }
        
        return true
    }
    
    /**
     *
     */
    protected static def filterTpedFiles(outBase){
        def test = "--recode ${filtOpts} --transpose "
        def comm = plinkTestTped(test, outBase+'.tped', outBase+'.tfam', outBase+'.filter')
        
        if( comm.execute().waitFor()!=0 ){
            System.err.println("Error recoding plink files [${outBase}]")
            return null
        }
        
        "mv ${outBase}.filter.tped ${outBase}.tped".execute().waitFor()
        "mv ${outBase}.filter.tfam ${outBase}.tfam".execute().waitFor()
        
        return true
    }
    
    /**
     * download a SNPs genotypes file from HapMap
     * 
     * @return the downloaded file path
     */
    protected def downloadSNPs(group){
        def file = hapmap_file.replace('{group}', group)
        file = file.replace('{chr}', chr)
        
        if( !Utils.download( HAPMAP_URL+file, workDir+file) ){
            System.err.println("Cannot download ${HAPMAP_URL+file}")
            return null
        }
        
        return workDir+file
    }
    
    /**
     *
     */
    protected def clearHapMapFiles(group){
        def file = hapmap_file.replace('{group}', group)
        file = file.replace('{chr}', chr)
        "rm -f ${workDir+file}".execute().waitFor()
        
        [trios_phased_file, duos_phased_file, unrel_phased_file].each{
            file = phasedFileName(it, group)
            "rm -f ${workDir+file}".execute().waitFor()
        }
    }
    
    /**
     *
     */
    protected def String phasedFileName(file, group){
        def lgroup = group.toLowerCase()
        file = file.replace('{lgroup}', lgroup)
        file = file.replace('{chr}', chr)
        return file
    }
    
    /**
     * download SNPs phasing files from HapMap
     * 
     * @return
     */
    protected def downloadPhasingSNPs(group){
        def urlBase = HAPMAP_PHASED_URL.replace('{group}', group)
        boolean error = false
        
        [(trios_phased_file):'TRIOS/', (duos_phased_file):'DUOS/', 
            (unrel_phased_file):'UNRELATED/'].each{ file, dir->
            file = phasedFileName(file, group)
            
            if( !Utils.download( urlBase+dir+file, workDir+file) ){
                System.err.println("Cannot download ${urlBase+dir+file}")
                error = true
            }
        }
        
        return !error
    }
    
    /**
     * parses the given line and extracts subjects ids
     * 
     * @return a map of {index,subject_id} values
     */
    private def getSubjectIdxs(String line, int start, int end, int step){
        def toks = line.split("\\s")
        def subjIdxs = [:]
        
        end = (end<1) ? toks.length : end
        
        start.step(end, step){ subjIdxs[it] = toks[it] }
        
        return subjIdxs
    }
    
    /**
     *
     */
    protected def parseHapMap(file, int start, int end){
        def reader = Utils.createReader(file)
        
        def header = reader.readLine()
        def subjIdxs = getSubjectIdxs(header, 11, 0, 1)
        
        genotypes.clear()
        subjIdxs.each{ idx, id-> genotypes[id] = new Genotype(subject: id) }
        
        reader.splitEachLine("\\s"){ toks->
            int pos = toks[3] as Integer
            
            if( pos>=start && pos<=end ){
                def currSnp = new SNPData(id:toks[0], chr:toks[2], 
                    alleles:toks[1][0]+toks[1][2], position:pos)
                
                hmSortSnps << currSnp
                hmSnps[toks[0]] = currSnp
                
                (11..toks.size()-1).each{
                    genotypes[subjIdxs[it]].snps[currSnp.id] = toks[it]
                }
            }
        }
        
        reader.close()
    }
    
    /**
     *
     */
    protected def parseHapMapPhasing(int start, int end, String group){
        [trios_phased_file, duos_phased_file, unrel_phased_file].each{ file->
            parseHapMapPhasing(workDir+phasedFileName(file, group), start, end)
        }
    }
    
    /**
     *
     */
    protected def parseHapMapPhasing(file, int start, int end){
        def reader = Utils.createReader(file)
        
        def header = reader.readLine()
        def subjIdxs = getSubjectIdxs(header, 2, 0, 2)
        
        subjIdxs.keySet().each{ idx->
            def id = subjIdxs[idx]
            id = id.substring(0, id.indexOf('_'))
            subjIdxs[idx] = id
            haplotypes[id] = new Genotype(subject: id) 
        }
        
        //line format: rsId position All1 All2 ...
        reader.splitEachLine("\\s"){ toks->
            int pos = toks[1] as Integer
            
            if( pos>=start && pos<=end ){
                if( !hmHaplo.containsKey(toks[0]) ){
                    def currSnp = new SNPData(id:toks[0], chr:chr, position:pos)
                    hmSortHaplo << currSnp
                    hmHaplo[toks[0]] = currSnp
                }
                
                2.step(toks.size(), 2){ 
                    haplotypes[subjIdxs[it]].snps[toks[0]] = toks[it]+toks[it+1] 
                }
            }
        }
        
        reader.close()
    }
    
    /**
     *
     * @param subjects : a list of subjects ids
     */
    static vcfToPed(vcfFile, tbiFile, String chr, int start, int end, 
        String outPref, boolean phased, subjects, tmpDir, plinkFormat = '--plink') {
        // add individuals to keep
        def keep = (subjects) ? subjects.sum{"--indv ${it} "} : ''
        chr = (chr.startsWith('chr')) ? chr.substring(3) : chr
        
        if( tbiFile ){
            //if tbiFile present slice vcf using tabix
            def params = [TABIX, '-h', vcfFile, "${chr}:${start}-${end}"]
            params = params.collect{it as String}
            
            if( vcfFile.startsWith('ftp://') )
                vcfFile = tmpDir + SLICE_VCF
            else
                vcfFile = vcfFile.substring(0, vcfFile.lastIndexOf('/')+1)+SLICE_VCF
            
            if( !procBuildRun(params, tmpDir, vcfFile) ){
                System.err.println("Cannot convert vcf to plink")
                return null
            }
            
            //compress new vcf file
            "gzip -f ${vcfFile}".execute().waitFor()
            vcfFile += '.gz'
        }
        
        def command = 
            "${VCFTOOLS} --gzvcf ${vcfFile} --out ${outPref} ${plinkFormat} "+
            "--chr ${chr} --from-bp ${start} --to-bp ${end} --remove-indels ${keep} "+
            "${(phased) ? '--phased' : ''}"
        
        if( command.execute().waitFor()!=0 ){
            System.err.println("Cannot convert vcf to plink")
            return null
        }
        
        return outPref
    }
    
    /**
     *
     */
    protected static vcfToTped(vcfFile, tbiFile, String chr, int start, int end, 
        String outPref, boolean phased, subjects, tmpDir) {
        
        return vcfToPed(vcfFile, tbiFile, chr, start, end, outPref, phased, 
            subjects, tmpDir, '--plink-tped')
    }
    
    /**
     * used in trans-eqtl calculations; converts a dir of .vcf files into .tped 
     * and .bgl also fills genotypes and imputed snps.
     */
    protected def vcfDirToTpedAndBgl(vcfDir, String outPref, subjects, snpsFile) {
        // add individuals to keep
        def keep = (subjects) ? subjects.sum{"--indv ${it} "} : ''
        
        def vcfFiles = new File(vcfDir).list({d, f-> f ==~ regexVcf } as FilenameFilter).toList()
        
        //parse snps file {chr, snp_id}
        def snpsMap = [:] //key=chr, value=list of snps
        new File(snpsFile).eachLine{ line->
            def toks = line.split("\\s")
            
            if( !snpsMap.containsKey(toks[0]) )
                snpsMap[toks[0]] = []
                
            snpsMap[toks[0]] << toks[1]
        }
        
        vcfFiles.each{ vcfFile ->
            def chr = (vcfFile =~ regexVcf)[0][1]
            
            if( snpsMap.containsKey(chr) ) {
                println "Found ${snpsMap[chr].size()} SNPs in ${chr}"
                def snps = snpsMap[chr].sum{"--snp ${it} "}
                
                chr = (chr.startsWith('chr')) ? chr.substring(3) : chr
                def outPrefChr = "${outPref}_chr${chr}"

                println "Processing SNPs of chr${chr}: ${vcfFile}"

                def command = "${VCFTOOLS} --gzvcf ${vcfDir+vcfFile} --out ${outPrefChr} --phased "+
                    "--plink-tped ${keep} --remove-indels --chr ${chr} ${snps}"

                if( command.execute().waitFor()!=0 ){
                    System.err.println("Cannot convert vcf chr${chr} to plink")
                    return null
                }

                if( filterTpedFiles(outPrefChr) ){
                    generateBgl(outPrefChr, true)
                }
                else{
                    println "Error filtering tped file ${outPrefChr}. Maybe no SNPs were selected."
                }
            }
        }
        
        return outPref
    }
        
    /**
     * run a process using ProcessBuilder
     */
    static boolean procBuildRun(params, workDir, outFile){
        ProcessBuilder pb = new ProcessBuilder(params)
        if( workDir )
            pb.directory(new File(workDir))

        //run command
        Process pro= pb.start()

        //get command output
        def istr = pro.getInputStream()
        def output = new FileOutputStream(outFile)
        
        output << istr
        
        try {
            //wait for process termination
            pro.waitFor();
        } catch (InterruptedException ex) {
        }
        
        istr.close()
        output.flush()
        output.close()
        
        if( pro.exitValue()!=0 )
            return false
        
        return true
    }
    
    /**
     * writes previously parsed hapmap genotypes to ped/map file
     * 
     * @param subjects : a list of subjects ids
     */
    protected def hapMapToPed(String outPref, genoMap, listSnps){
        //write ped file
        def writer = new PrintWriter(outPref+'.ped')
        
        genoMap.each{ key, geno->
            writer.print("${key} ${key} 0 0 0 0")//header   
            listSnps.each{ //genotype
                def alleles = geno.snps[it.id]
                writer.print(" ${alleles[0]!='N' ? alleles[0] : '0'} ${alleles[1]!='N' ? alleles[1] : '0'}")
            }
            writer.print("\n")
        }
        
        writer.close()
        
        //write map file
        writer = new PrintWriter(outPref+'.map')
        listSnps.each{
            writer.println("${it.chrNum} ${it.id} 0 ${it.position}")
        }
        writer.close()
    }
    
    /**
     * builds a plink command call
     * 
     * @param plink : plink executable
     * @param test : --assoc, --model, ...
     */
    public static String plinkTest(test, ped, map, output){
        return "${PLINK} --noweb --ped ${ped} --map ${map} ${test} --out ${output}"
    }
    
    public static String plinkTestTped(test, tped, tfam, output){
        return "${PLINK} --noweb --tped ${tped} --tfam ${tfam} ${test} --out ${output}"
    }
    
    /**
     * generates beagle genotype file associated to tped/tfam file and populates
     * imputedSnps map
     */
    private generateBgl(outPref, boolean addMode = false){
        def writer = new PrintWriter(outPref+'.bgl')
        
        //write header and create genotypes
        if( !addMode )
            genotypes.clear()
        
        def subjIdxs = [:]
        int count = 0
        
        writer.print("I\tid")
        def reader = Utils.createReader( new File(outPref+'.tfam') )
        reader.splitEachLine("\\s"){ toks->
            writer.print("\t${toks[1]}\t${toks[1]}")
            
            if( !genotypes[toks[1]] )
                genotypes[toks[1]] = new Genotype(subject: toks[1])
                
            subjIdxs[count++] = toks[1]
        }
        writer.print("\n")
        reader.close()
        
        //write genotypes
        reader = Utils.createReader( new File(outPref+'.tped') )
        
        reader.eachLine{ line->
            def toks = line.split("\\s",5)
            def snpId = toks[1]
            writer.println("M\t${snpId}\t${toks[4]}")
            
            //read genotypes
            (0..subjIdxs.size()-1).each{
                int idx = it*4
                genotypes[subjIdxs[it]].snps[snpId] = (toks[4][idx]+toks[4][idx+2])
            }
            
            //get snp alleles
            def alleles = [] as Set
            for(int i=0; i<toks[4].length(); i++){ 
                if( toks[4][i] in ['A','T','G','C'] )
                    alleles << toks[4][i]
                    
                if(alleles.size()==2)
                    break
            }
            
            if( alleles.size()<2 ) {
                println "Warning: snp ${snpId} chr ${toks[0]} is homozygotic in all subjects {${alleles.sum()}}"
            }
            
            //add snp to imputed snps and imputedSort list
            imputedSnps[snpId] = new SNPData( id:snpId, chr:toks[0], 
                position:toks[3] as Integer, alleles:alleles.sort().sum() )
            imputedSort << imputedSnps[snpId]
        }
        
        reader.close()
        writer.close()
    }
    
    /**
     * generates beagle reference panel file associated to tped/tfam file,
     * also generates markers file
     */
    private generatePhasedBgl(outPref){
        def writer = new PrintWriter(outPref+'.bgl')
        
        //write header
        writer.print("I\tid")
        def reader = Utils.createReader(new File(outPref+'.tfam'))
        reader.splitEachLine("\\s"){ toks->
            writer.print("\t${toks[1]}\t${toks[1]}")
        }
        writer.print("\n")
        reader.close()
        
        //write genotypes
        reader = Utils.createReader(new File(outPref+'.tped'))
        def writerMark = new PrintWriter(outPref+'.mark')
        
        def tpedLines = [:] as TreeMap
        // read .tped lines
        reader.eachLine{ line->
            def toks = line.split("\\s",5)
            tpedLines[toks[1]] = toks[4]
        }
        
        imputedSort.each{ snpdata->
            if( tpedLines.containsKey(snpdata.id) ){
                //write bgl line
                writer.println("M\t${snpdata.id}\t${tpedLines[snpdata.id]}")
                //write marker line (snpid,position,allele1,allele2)
                writerMark.println("${snpdata.id}\t${snpdata.position}\t${snpdata.alleles[0]}\t${snpdata.alleles[1]}")
            }
            else
                println "generatePhasedBgl: SNP ${snpdata.id} not present in merged genotypes set"
        }
        
        reader.close()
        writer.close()
        writerMark.close()
    }
    
    /**
     * 
     */
    private generateMergedBgl(outPref, tpedPref){
        def writer = new PrintWriter(outPref+'.bgl')
        
        //write header
        def subjIdxs = [:]
        int count = 0
        
        writer.print("I\tid")
        def reader = Utils.createReader(new File(tpedPref+'.tfam'))
        reader.splitEachLine("\\s"){ toks->
            writer.print("\t${toks[1]}\t${toks[1]}")
            subjIdxs[count++] = toks[1]
        }
        writer.print("\n")
        reader.close()
        
        //write genotypes
        reader = Utils.createReader(new File(tpedPref+'.tped'))
        
        reader.eachLine{ line->
            def toks = line.split("\\s",5)
            def snpId = toks[1]
            writer.print("M\t${snpId}")
            
            //write genotypes
            (0..subjIdxs.size()-1).each{
                def alleles = genotypes[subjIdxs[it]].snps[snpId]
                writer.print("\t${alleles[0]}\t${alleles[1]}")
            }
            
            writer.print("\n")
        }
        
        reader.close()
        writer.close()
    }
    
    /**
     * parses bgl.phased file generated by beagle
     */
    protected def parseBgl(fileBgl, double r2Cut){
        //read snps
        def fileDose = fileBgl.replace('.phased','.dose')
        def reader = Utils.createReader(new File(fileDose))
        
        reader.readLine()
        reader.eachLine{ line->
            def toks = line.split("\\s",4)
            // set alleles
            imputedSnps[toks[0]].alleles = [toks[1],toks[2]].sort().sum()
        }
        reader.close()
        
        //discard snps by r2 cutoff
        def fileR2 = fileBgl.replace('.phased.gz','.r2')
        reader = Utils.createReader(new File(fileR2))
        
        
        reader.splitEachLine("\\s"){ toks->
            try{ 
                double r2 = toks[1] as Double
                if( r2 < r2Cut )
                    discardedSnps << toks[0]
            } catch(e){
                discardedSnps << toks[0]
            }
        }
        reader.close()
        
        println "${discardedSnps.size()} snps of ${imputedSnps.size()} filtered by R2 cutoff of ${r2Cut}"
        
        //read imputed genotypes
        reader = Utils.createReader(new File(fileBgl))
        
        def header = reader.readLine()
        def subjIdxs = getSubjectIdxs(header, 2, 0, 2)
        
        reader.splitEachLine("\\s"){ toks->
            def currSnp = imputedSnps[toks[1]]
            
            if( !(currSnp.id in discardedSnps) ){
                2.step(toks.size(), 2){
                    genotypes[subjIdxs[it]].snps[currSnp.id] = toks[it]+toks[it+1]
                }
            }
        }
        
        reader.close()
    }
    
    /**
     *
     */
    protected def writeGenotypes(file){
        def writer = new PrintWriter(file)
        def subjects = genotypes.keySet().sort()
        
        //write header
        writer.print("SNP")
        subjects.each{ writer.print("\t${it}") }
        writer.print("\n")
        
        imputedSnps.keySet().each{
            writer.print("${it}")
            subjects.each{ subj->
                def alleles = genotypes[subj].snps[it]
                alleles = [alleles[0],alleles[1]].sort().sum()
                writer.print("\t${alleles}")
            }
            writer.print("\n")
        }
        
        writer.close()
    }
    
    /**
     * load and run beagle
     */
    protected static def runBeagle(beagleJar, args){
//        URL jar = SNPManager.getClassLoader().getSystemResource(BEAGLE_JAR)
//        URLClassLoader ucl = new URLClassLoader([jar] as URL[]);
//        def phaseclass = ucl.loadClass("phaser.PhaseMain");
//        phaseclass.main(args as String[])
        String javahome = System.getProperty("java.home")+File.separator+"bin"+File.separator
        def comm = "${javahome}java -Xmx1g -jar ${beagleJar} ${args.sum{it+' '}}"
        println "Running: ${comm}"
        def proc = comm.execute()
        proc.consumeProcessOutputStream(System.out)
        proc.consumeProcessErrorStream(System.err)
        proc.waitFor()
    }
    
    protected static String installBeagle(outdir){
        def istr = SNPManager.getClassLoader().getSystemResourceAsStream(BEAGLE_JAR)
        def outstr = new FileOutputStream(outdir+BEAGLE_JAR)
        StorageManager.writeInputToOutput( istr, outstr )
        return outdir+BEAGLE_JAR
    }
    
    /**
     *
     */
    static downloadFromS3(s3url, outdir){
        def bucket = StorageManager.getBucket(s3url)
        def s3object = StorageManager.removePrefix(s3url)
        def file = outdir + s3object.substring(s3object.lastIndexOf('/')+1)
        
        if( S3Manager.downloadFile( bucket, s3object, file)==null )
            throw new StorageException("Error downloading file: ${s3object}\n bucket: ${bucket}, dst: ${file}")
 
        return file
    }
    
}

