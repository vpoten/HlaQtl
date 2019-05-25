/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table

import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.IntColumn
import tech.tablesaw.api.DoubleColumn
import tech.tablesaw.io.csv.CsvWriteOptions
import groovy.cli.picocli.CliBuilder
import groovy.cli.Option
import groovy.cli.Unparsed

import org.msgenetics.hlaqtl.eqtl.LDCalculator
import org.msgenetics.hlaqtl.eqtl.SNPManager
import org.msgenetics.hlaqtl.eqtl.GTExEqtl


/**
 *
 * @author victor
 */
class GTExSearcher {
    
    static final String ENSEMBL_REST_API = 'http://grch37.rest.ensembl.org'
    
    /** list of query snps (rs ids) */
    List<String> queryIds = []
    
    /** Size of region around query SNPs */
    int snpRegionSize = 10000000
    
    /** GTEx eqtl data */
    String gtexDir = null
    
    /** Working directory */
    String workDir = null
    
    /** 1000 genomes (vcf + tbi) directory */
    String genomesDir = null
    
    /** GTEx tissues where to filter the eqtls (all by default) */
    List<String> tissues = null
    
    /** eqtl p-value threshold to filter best eqtls */
    double eqtlThr = 0.05d
    
    /** LD threshold to filter results */
    double ldThr = 0.5d
    
    /** Keep cache of tped files */
    boolean useCache = true
    
    /** List of subjects to use */
    List<String> subjects
    
    /** Number of threads to use in LD calculation */
    int nThreads = 4
    
    /** Allowed chromosomes */
    Set<String> chrAllowed = null;
    
    /**
     * Main method, parse args and perform the calculation
     */ 
    public static void main(args) {
        def instance =  createFromArgs(Arrays.copyOfRange(args, 1, args.length))
        if (instance != null) {
            instance.perform()
        }
    }
    
    static CliBuilder cliBuilder() {
        def cli = new CliBuilder(name: Main.COMM_GTEX_SEARCH)
        cli._(longOpt: 'help', args: 0, 'display help')
        cli._(longOpt: 'genomesDir', argName:'path', args: 1, 'directory where 1000 genomes vcf files resides', required: true)
        cli._(longOpt: 'gtexDir', argName:'path', args: 1, 'directory where eQTL GTEx results resides', required: true)
        cli._(longOpt: 'snps', argName:'path', args: 1, 'file containing the rs ids of the query SNPs', required: true)
        cli._(longOpt: 'workDir', argName:'path', args: 1, 'working directory', required: true)
        cli._(longOpt: 'populations', argName:'code', valueSeparator:',', args: '+', defaultValue: 'CEU', 'populations where to get the subjects')
        cli._(longOpt: 'eqtlThr', argName:'thr', args: 1, defaultValue: '0.05', 'eQTL p-value threshold', type: Double)
        cli._(longOpt: 'ldThr', argName:'thr', args: 1, defaultValue: '0.5', 'LD result thresholdd', type: Double)
        cli._(longOpt: 'regionSize', argName:'len', args: 1, 'SNP region size', defaultValue: '10000000', type: Integer)
        return cli
    }
    
    /**
     * Create instance from commandline args
     */ 
    static GTExSearcher createFromArgs(args) {
        def cli = cliBuilder()                  
        def options = cli.parse(args)
        
        if (options == null) {
            return null
        }
        
        def instance = new GTExSearcher()
        instance.setEqtlThr(options.eqtlThr)
        instance.setLdThr(options.ldThr)
        instance.setSnpRegionSize(options.regionSize)
        instance.loadQuerySnpsFromFile(options.snps)
        instance.setSubjects(options.populations)
        return instance
    }
    
    /**
     *
     */
    def perform() {
        // Load snp info from ensembl rest api
        EnsemblRestClient client = new EnsemblRestClient(ENSEMBL_REST_API, 15, 200)
        List<SNPData> snpQuery = client.getSnps(queryIds, 'human')
        
        // apply allowed chromosomes restriction
        if (chrAllowed != null) {
            snpQuery = snpQuery.findAll{it.chr in chrAllowed}
        }
        
        // Load best eqtls from GTEx data
        Table bestEqtls = GTExEqtl.getBestEqtlsAllTissues(gtexDir, eqtlThr);
        
        // Build regions around snps in query list
        def chrRegions = [:] as TreeMap
        
        snpQuery.each { snp ->
                if ( !(snp.chrNum in chrRegions) ) {
                    chrRegions[snp.chrNum] = []
                }
                def start = (snp.position - snpRegionSize) < 1 ? 1 : snp.position - snpRegionSize
                chrRegions[snp.chrNum] << ['start': start, 'end': snp.position + snpRegionSize, 'snp': snp]
            }
        
        // sort chr regions by start position
        for(String chr : chrRegions.keySet()) {
            Collections.sort(chrRegions[chr], [compare: {a, b -> a.start <=> b.start}] as Comparator)
        }
        
        // Associate eqtls with regions
        chrRegions.each { chr, regions ->
            regions.each { region ->
                Table result = GTExEqtl.filterByRegion(bestEqtls, chr, region.start, region.end)
                region['eqtls'] = result
            }
        }
        
        if(useCache == false || !checkExistsTped(chrRegions)) {
            // Obtain tped files from 1000genomes vcfs
            extractTpedFromVcfs(chrRegions)
        }
        
        //  Load SNPData from tped files for all snps: query + eqtl
        def snpsData = [:] as TreeMap
        
        chrRegions.each { chr, regions ->
            def tpedFile = new File(buildChrDir(chr), "${SNPManager._1000G_PED}.tped").absolutePath
            def snps = [] as TreeSet
            
            regions.each { region->
                // add the snp that leads the region and all the snps in associated eqtls
                snps << region.snp.id
                Table eqtls = region.eqtls
                def regionSnps = ((StringColumn) eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13')).
                    asSet().
                    findAll{ it!=null && it.startsWith('rs') }
                region['region_snps'] = regionSnps
                snps += regionSnps
            }
            
            println "Load genotypes from tped file chr${chr}: ${new Date()}"
            SNPData.createFromTped(tpedFile, snps, snpsData)
        }
        
        // write regions report
        reportRegionsStats(chrRegions, snpsData)
        
        // LD calculation between query snps and eqtl snps in region
        println "\nCalculating LD in ${nThreads} threads: ${new Date()}\n"
        regionLDCalc(chrRegions, snpsData)
        
        // Final result
        Table finalTableFiltered = Table.create("GTEx_eqtl_search_filter_LD")
        Table finalTableAll = Table.create("GTEx_eqtl_search_all_LD")

        chrRegions.each { chr, regions ->
            regions.each { region->      
                if (!('ld_results' in region)) {
                    region['ld_results'] = [:]
                }
                
                // get snpIds with rSq > threshold
                def ldResultsPass = region['ld_results'].findAll({ snpId, rSq -> rSq > ldThr})
                
                // create the final result table for the region: eqtl table + extra columns:
                // extra columns: region_snp, region_snp_pos, ld_rsq
                
                StringColumn regionSnps = (StringColumn) region.eqtls.stringColumn('rs_id_dbSNP147_GRCh37p13')
                Table eqtlsSubset = region.eqtls.where(regionSnps.isIn(ldResultsPass.keySet()))
                
                finalTableFiltered = appendToResultTable(region, eqtlsSubset, finalTableFiltered)
                finalTableAll = appendToResultTable(region, region.eqtls.copy(), finalTableAll)
            }
        }
        
        [finalTableAll, finalTableFiltered].each { it ->
            if (it.rowCount() > 0) {
                writeResultTable(it)
            }
            else {
                println "Empty result: No table ${it.name()} will be written to disk.\n"
            }
        }
    }
    
    /**
     * 
     */ 
    private Table appendToResultTable(region, Table eqtlsSubset, Table destination) {
        if (eqtlsSubset.rowCount() == 0) {
            return destination
        }
        
        StringColumn resultSnps = (StringColumn) eqtlsSubset.stringColumn('rs_id_dbSNP147_GRCh37p13')
        // create extra columns
        def snpCol = StringColumn.create('region_snp', (0..eqtlsSubset.rowCount()-1).collect{region.snp.id})
        def snpPosCol = IntColumn.create('region_snp_pos', (0..eqtlsSubset.rowCount()-1).collect{region.snp.position} as int [])
        def ldRsqCol = DoubleColumn.create('ld_rsq', resultSnps.asList().collect{region['ld_results'][it] ?: -1})

        // add the columns to the result table
        eqtlsSubset = eqtlsSubset.addColumns(snpCol, snpPosCol, ldRsqCol)

        // append to final table
        if (destination.columnCount() == 0) {
            destination = destination.addColumns(eqtlsSubset.columnArray())
        }
        else {
            destination = destination.append(eqtlsSubset)
        }
        
        return destination
    }
    
    /**
     *
     */
    private def regionLDCalc(chrRegions, snpsData) {
        // collect all regions
        def allRegions = chrRegions.collect({chr, regions -> regions}).flatten()
        
        // closure for LD calculation (threaded)
        def lDThreadCalc = { thread, totalThrs ->
            // assign index to region to divide jobs among threads
            allRegions.eachWithIndex{ region, i-> 
                if ((i % totalThrs) == thread) {
                    region['ld_results'] = [:] as TreeMap
                    def current = region.snp.id
                    def data1 = snpsData[current]

                    if (data1 != null) {
                        region['region_snps'].each { snpId ->
                            if( current != snpId ) {
                                def data2 = snpsData[snpId]
                                if (data2 != null) {
                                    Double rSq = LDCalculator.calcLD(data1, data2)
                                    if( rSq!=null ) {
                                        region['ld_results'][snpId] = rSq
                                    }
                                }
                                else {
                                    System.err.println("SNP ${snpId} in region lead by ${current} (${region.snp.locus}) has no data")
                                }
                            }
                        }
                    }
                    else {
                        System.err.println("Region SNP ${current} (${region.snp.locus}) has no data")
                    }
                }
            }
        }// end lDThreadCalc closure
        
        // run lDThreadCalc in threads
        nThreads = 4
        
        def closures = [
            {lDThreadCalc(0, nThreads)},
            {lDThreadCalc(1, nThreads)},
            {lDThreadCalc(2, nThreads)},
            {lDThreadCalc(3, nThreads)}
        ]
        
        Utils.runClosures(closures, nThreads )
    }
    
    /**
     * Set subjects to use in calculations
     */
    def setSubjects(populations) {
        subjects = IGSRSamples.create().getSubjects(populations)
    }
    
    /**
     * Load query snp ids from text file
     */
    def loadQuerySnpsFromFile(file) {
        def reader = new File(file).newReader()
        queryIds.clear()
        
        reader.eachLine { line ->
            queryIds << line.trim()
        }
        
        reader.close()
    }
    
    /**
     * Build a directory for a chromosome within workDir and returns the created path
     */
    def buildChrDir(String chr) {
        if( !chr.startsWith('chr') ){
            chr = 'chr' + chr
        }
        def chrDir = new File(workDir, chr)
        String chrPath = chrDir.absolutePath
        
        if (!chrDir.exists()) {
            Utils.createDir(chrPath)
        }
        
        return chrPath
    }
    
    /**
     * Set allowed chromosomes to use in calculations
     */
    def setChrAllowed(chrs) {
        chrAllowed = chrs as Set
    }
    
    /**
     * Write csv result table to disk
     */
    private void writeResultTable(Table result) {
        def builder = CsvWriteOptions.builder(new File(workDir, "${result.name()}_${new Date().format('yyyyMMddHHmmss')}.csv")).
            header(true).
            separator((char)'\t')
        
        CsvWriteOptions options = builder.build()
        result.write().csv(options)
    }
    
    /**
     * Checks that tped files exists in working directory
     */
    private boolean checkExistsTped(chrRegions) {
        return chrRegions.keySet().collect{new File(buildChrDir(it), "${SNPManager._1000G_PED}.tped")}.every{it.exists()}
    }
    
    /**
     * Extract tped from vcfs in threads
     */
    private void extractTpedFromVcfs(chrRegions) {
        
        // closure for tped extraction
        def tpedFromVcfThreadCalc = { thread, totalThrs ->
            chrRegions.keySet().eachWithIndex { chr, i ->
                if ((i % totalThrs) == thread) {
                    def regions = chrRegions[chr]
                    int start = regions.min{ it.start }.start
                    int end = regions.max{ it.end }.end
                    def locusStr = "${chr}:${start}-${end}"
                    def vcfTpl = (chr == 'X') ? SNPManager.S3_VCF_FILE_X : SNPManager.S3_VCF_FILE
                    def vcfFile = new File(genomesDir , vcfTpl.replace('{chr}', "chr${chr}")).absolutePath
                    def chrDir = buildChrDir(chr) + '/'
                    // generate tped files in a separate directory for each chromosome
                    println "Load snp data from vcf file chr${chr}: ${new Date()}"
                    SNPManager.loadSNPData(subjects, chrDir, vcfFile, [], locusStr, true, null)
                }
            }
        }
        
        // run tpedFromVcfThreadCalc in threads
        def closures = [
            {tpedFromVcfThreadCalc(0, 2)},
            {tpedFromVcfThreadCalc(1, 2)}
        ]
        
        Utils.runClosures(closures, 2)
    }
    
    /**
     * Writes a report with region info to disk
     */
    private void reportRegionsStats(chrRegions, snpsData) {
        def writer = new File(workDir, 'regions_summary.tsv').newPrintWriter()
        writer.println('chr\tposition\tsnp\treg_start\treg_end\tnum_eqtl\tnum_snps\tmissing_genotypes\trate_missing')
        
        chrRegions.each { chr, regions ->
            regions.each { region->
                region['missing_snps'] = region['region_snps'].findAll{!(it in snpsData)}
               
                writer.print("${chr}\t${region.snp.position}\t${region.snp.id}\t")
                writer.print("${region.start}\t${region.end}\t")
                writer.print("${region.eqtls.rowCount()}\t${region['region_snps'].size()}\t${region['missing_snps'].size()}\t")
                writer.println("${region['missing_snps'].size()/(double)region['region_snps'].size()}")
            }
        }
        
        writer.close()
    }
}

