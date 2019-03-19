/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.msgenetics.hlaqtl.eqtl.FisherTestResult
import org.ngsutils.Utils
import org.ngsutils.variation.SNPData
import org.ngsutils.annotation.genome.AnnotationTrack as AnT
import org.ngsutils.annotation.genome.AnnotationFactory as AnF
import org.ngsutils.annotation.genome.FeatureIndex
import org.ngsutils.annotation.genome.regulation.TFBMotif
import org.ngsutils.stats.FisherExact
import net.sourceforge.jFuzzyLogic.FIS
import org.apache.commons.math3.stat.inference.TestUtils

/**
 * 
 */
class RegulationPattern {
    SNPData snp
    String line
    Double score = null
    
    def cellMap = null //key=cell_line, value=map(key=histone, value=count)
    def regulationMap = null //key=histone, value=[#inCell,#otherLines,sumOther]
    
    def tfsInCl = [] as TreeSet //TFBS list from Uniform TFs database
    def motifs = [:] as TreeMap //key=motif, value=position_weigth
    def transFac = [] as TreeSet //TFBS list from clustered TFs database
    
    static final MISSING = ['null','NA','-']
    
    protected def createRegTrkMap = {
        def map = [:]
        AnF.regulationTracks.each{ map[it] = 0 }
        return map
    }
    
    // create cell line map for regulation tracks
    protected def createCellRegMap = {
        def map = [:]
        AnF.cellLines.each{ map[it] = createRegTrkMap() }
        return map
    }
    
    
    /**
     * create cell line count map for histones and TFs
     */
    def createCellMap() {
        cellMap = createCellRegMap()
    }
    
    /**
     * 
     */
    def update(String _line, String featName) {
        if( featName in AnF.regulationTracks ){
            cellMap[_line][featName]++ //is a DNase or histone
        }
        else{
            if( line==_line ){
                tfsInCl << featName
            }
        }
    }
    
    
    /**
     * 
     */ 
    static tabHeader(writer) {
        writer.print("snp\tlocus\tline\t")
        AnF.regulationTracks.each{ writer.print("${it}_in_CL\t${it}_in_other\t${it}_in_other_sum\t") }
        writer.println("TFs_in_CL\tTFs\tMotifs")
    }
    
    /**
     * 
     */ 
    def printTab(writer) {
        writer.print("${snp.id}\t${snp.locus}-${snp.position+1}\t${line}\t")
        AnF.regulationTracks.each{ trk->
            if( cellMap ) {
                int sumOther = AnF.cellLines.sum{ (it==line) ? 0 : cellMap[it][trk] }
                int sumLines = AnF.cellLines.sum{ (it==line) ? 0 : ((cellMap[it][trk]>0) ? 1 : 0) }
                writer.print("${cellMap[line][trk]}\t${sumLines}\t${sumOther}\t")
            }
            else{
                writer.print("${regulationMap[trk][0]}\t${regulationMap[trk][1]}\t${regulationMap[trk][2]}\t")
            }
        }
        
        if( motifs && !transFac){ System.println("Warning: in pattern ${snp.id}; motifs and no TFs.") }
        
        writer.print( "${tfsInCl?.sum{it+','}}\t" )
        writer.print( "${transFac.sum{it+','}}\t" )
        writer.println( motifs.keySet().sum{"${it}=${motifs[it]},"} )
    }
    
    /**
     * 
     */
    Integer getValue(trk, boolean inOther) {
        if( regulationMap ) {
            return regulationMap[trk][(inOther ? 1 : 0)]
        }
        else if( cellMap ) {
            if( inOther ) {
                return AnF.cellLines.sum{ (it==line) ? 0 : ((cellMap[it][trk]>0) ? 1 : 0) }
            }
            else {
                return cellMap[line][trk]
            }
        }
        
        return null
    }
    
    /**
     * @return a map with key=snpId and value=RegulationPattern
     */ 
    static Map loadTab(file, set=null) {
        def createRegMap = {
            def map = [:]
            AnF.regulationTracks.each{ map[it] = [0,0,0] }
            return map
        }
        
        
        def reader = Utils.createReader(new File(file))
        reader.readLine() //skip header
        def map = [:] as TreeMap
        
        reader.splitEachLine("\t"){ toks->
            if( set==null || (toks[0] in set) ) {
                def mat = (toks[1]=~Utils.locusRegex)[0]

                def snp = new SNPData(id:toks[0], chr:mat[1], position:(mat[2] as Integer))
                RegulationPattern patt = new RegulationPattern(snp:snp, line:toks[2])
                patt.regulationMap = createRegMap()

                AnF.regulationTracks.eachWithIndex{ trk, i->
                    int base = 3 + 3*i
                    (0..2).each{ patt.regulationMap[trk][it] = toks[base+it] as Integer }
                }

                int tfIdx = 3 + 3*AnF.regulationTracks.size()

                //tf in cell  line
                if( !(toks[tfIdx] in MISSING) ) {
                    patt.tfsInCl += (toks[tfIdx].split(',') as List)
                }
                
                //tf clustered
                if( !(toks[tfIdx+1] in MISSING) ) {
                    patt.transFac += (toks[tfIdx+1].split(',') as List)
                }

                //motifs
                if( !(toks[tfIdx+2] in MISSING) ) {
                    (toks[tfIdx+2].split(',') as List).each{ tok->
                        int pos = tok.indexOf('=')
                        patt.motifs[tok.substring(0,pos)] = tok.substring(pos+1) as Double
                    }
                }

                map[toks[0]] = patt
            }
        }
        
        reader.close()
        return map
    }
    
}



/**
 *
 * @author victor
 */
class RegulationSnps {
	
    FIS fInference //fuzzy inference
    def normalization = [:]
    
    static final String OUT_GRP_RANK = 'group_rank.txt'
    static final String ANOVA_TEST_OUT = 'motifs_anova.out.txt'
    /**
     *
     * @param workDir
     * @param snpsInfoFile
     * @param cellLine = cell line where the assay is done
     */
    static List generatePatterns(String workDir, String snpsInfoFile, String cellLine) {
        // parse snps log file
        def snps = []
        def reader = new File(snpsInfoFile).newReader()
        reader.readLine() //skip header
        
        reader.splitEachLine("\\s"){ toks->
            def snp = new SNPData( //zero based genomic position (position-1)
                id:toks[0], chr:(toks[1]=='23' ? 'X' : toks[1]), position:((toks[2] as Integer)-1),
                minor:toks[4], major:toks[5], maf:toks[6] as Double )
            
            snps << snp
        }
        reader.close()
        generatePatterns(workDir, snps, cellLine)
    }
    
    /**
     *
     * @param workDir
     * @param snps = list of SNPData
     * @param cellLine = cell line where the assay is done
     */
    static List generatePatterns(String workDir, List snps, String cellLine) {
        def taxId = '9606'
        def assembly = 'hg19'
        
        // build histone + DNase tracks for each cell line
        def cellLinesTrks = AnF.broadHistoneAndUniDNaseTracks(assembly)
        def cellLinesIdxs = [:]
        
        // Load histone indexes
        cellLinesTrks.each{ name, tracks->
            FeatureIndex regIndex = 
                AnF.createIndex(workDir, taxId, tracks)
            cellLinesIdxs[name] = regIndex
        }
        
        // Load TFBinding indexes
        FeatureIndex tfIndex = AnF.createIndex(workDir, 
            taxId, AnF.ANN_UCSC_HG19_TF)
        FeatureIndex motifIndex = AnF.createIndex(workDir, 
            taxId, AnF.ANN_UCSC_HG19_Motif)
        
        // Load TFBinding Motifs
        def tfbMotifMap = TFBMotif.load(workDir, assembly, taxId)
        
        def patterns = []
        
        snps.each{ snp->
            String locus = "${snp.locus}-${snp.position+1}"//zero based
            RegulationPattern patt = new RegulationPattern(snp:snp, line:cellLine)
            patt.createCellMap()
            
            cellLinesIdxs.each{ line, regIndex->
                def feats = regIndex.getFeatsByPos(locus)
                feats?.each{ feat-> patt.update(line,feat.name) }
            }
            
            def feats = tfIndex.getFeatsByPos(locus, 10000)
            def motifs = motifIndex.getFeatsByPos(locus, 500)
            
            motifs?.each{ motif->
                Double prob = 
                    tfbMotifMap[motif.name]?.getMaxProb( snp.position-motif.location.start() )
                patt.motifs[motif.name] = (prob==null) ? 0.0d : prob
            }
            
            feats?.each{ feat->
                patt.transFac << feat.name
            }
            
            patterns << patt
        }
        
        return patterns
    }
    
    /**
     * writes regulation snp patterns to file
     */
    static writePatterns(file, patterns) {
        def writer = new PrintWriter(file)
        RegulationPattern.tabHeader(writer)
        patterns.each{ it.printTab(writer) }
        writer.close()
    }
    
    /**
     * Loads a fuzzy inference (.fcl) file
     */
    def loadInference(file=null) {
        def istr = (file) ? new File(file).newInputStream() : 
            RegulationSnps.getClassLoader().getSystemResourceAsStream("regulation.fcl")
        
        fInference = FIS.load(istr,true)
    }
    
    /**
     * evaluates pattern and sets its score
     */
    def setFuzzyScore(pattern) {
        def functionBlck = "regulation_onlyCL"
        
        AnF.regulationTracks.each{ trk->
            fInference.setVariable(functionBlck, trk, (pattern.getValue(trk,false)>0) ? 1 : 0)
        }
        
        if( functionBlck=='regulation') {
            AnF.regulationTracks.each{ trk->
                fInference.setVariable(functionBlck, trk+'_other', pattern.getValue(trk,true)*normalization[trk])
            }
        
            // clustered transFac
            fInference.setVariable(functionBlck, 'transFac', pattern.transFac.size()*normalization['transFac'])
        }
        else {
            //transFac in cell line
            fInference.setVariable(functionBlck, 'tfsInCl', pattern.tfsInCl.size()*normalization['tfsInCl'])
        }
        
        //motif
        fInference.setVariable(functionBlck, 'motif', (pattern.motifs) ? pattern.motifs.values().max() : 0)
        
        fInference.evaluate()
        
        
        pattern.score = fInference.getFunctionBlock(functionBlck).getVariable('score').getValue()
    }
    
    /**
     * Calculates normalization factor for fuzzy variables. This method must be
     * called before pattern evaluation
     */
    def calcNormalization( patterns ) {
        AnF.regulationTracks.each{ trk->
            normalization[trk] = 1.0d/(patterns.max{it.getValue(trk,true)}.getValue(trk,true))
        }
        
        normalization['tfsInCl'] = 1.0d/(patterns.max{it.tfsInCl.size()}.tfsInCl.size())
        normalization['transFac'] = 1.0d/(patterns.max{it.transFac.size()}.transFac.size())
    }
    
    /**
     *
     */
    static performGroupRanking(args) {
        String patternsFile = Main.readOption(args, Main.OPT_SNPS)
        String groupsFile = Main.readOption(args, Main.OPT_GROUP)
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        outdir = dirs[0]
        
        //load groups file
        def groups = [:] as TreeMap //key=group, value=list of snpIds
        def snps = [] as TreeSet
        
        def reader = Utils.createReader(new File(groupsFile))
        reader.readLine() //skip header
        reader.splitEachLine("\\s"){ toks->
            def group = groups[toks[0]]
            
            if( group==null ){
                group = []
                groups[toks[0]] = group
            }
            
            group << toks[1]
            snps << toks[1]
        }
        
        reader.close()
        
        // load patterns
        def patterns = RegulationPattern.loadTab(patternsFile)
        
        // do fuzzy inference
        def instance = new RegulationSnps()
        instance.loadInference()
        instance.calcNormalization( patterns.values() )
        
        snps.each{ instance.setFuzzyScore(patterns[it]) }
        
        //write output
        def scoreComparator = [
                compare:{a,b-> patterns[b].score<=>patterns[a].score }
            ] as Comparator
        
        def writer = new PrintWriter(outdir+OUT_GRP_RANK)
        writer.print("group\tscore\t")
        RegulationPattern.tabHeader(writer)
        
        groups.each{ grp, list->
            Collections.sort(list, scoreComparator)
            
            list.each{ 
                writer.print("${grp}\t${patterns[it].score}\t")
                patterns[it].printTab(writer)
            }
        }
        
        writer.close()
    }
    
    /**
     *
     */
    static performTransfacEnrichment(args) {
        String patternsFile = Main.readOption(args, Main.OPT_SNPS)
        String snpsFile = Main.readOption(args, Main.OPT_VALUE)
        String outdir = args.find{ it.startsWith(Main.OPT_OUTPUT_DIR) }
        
        def dirs = [outdir]
        dirs = Main.checkAndCleanDirOpts(dirs)
        outdir = dirs[0]
        
        //load snps
        def selectedSnps = [] as TreeSet
        def reader = Utils.createReader(new File(snpsFile))
        reader.splitEachLine("\\t"){ toks-> selectedSnps << toks[0] }
        reader.close()
        
        // load patterns
        def patterns = RegulationPattern.loadTab(patternsFile)
        
        // count TFs in all and selected distributions
        def transFacCounts = ['allTFs':[0,0]] as TreeMap
        def motifDists = ['allMotifs':[[],[],0,0]] as TreeMap
        
        // create counts and motif maps
        patterns.each{ snp, patt->
            patt.transFac.each{ 
                if(!transFacCounts.containsKey(it)) {
                    transFacCounts[it] = [0,0]
                }
            }
            
            patt.motifs.each{
                if(!motifDists.containsKey(it.key)) {
                    motifDists[it.key] = [[],[],0,0]
                }
            }
        }
        
        // count all_snps distribution
        patterns.each{ snp, patt->
            if( !(snp in selectedSnps) ){
                patt.transFac.each{  transFacCounts[it][0]++ }
                if( patt.transFac ){ transFacCounts['allTFs'][0]++ }
                
                patt.motifs.each{
                    motifDists[it.key][0] << it.value
                    motifDists[it.key][2]++ 
                }
                if( patt.motifs ){ 
                    motifDists['allMotifs'][0] += patt.motifs.values()
                    motifDists['allMotifs'][2]++
                }
            }
        }
        
        //count selected_snps distribution
        selectedSnps.each{ snpId->
            def patt = patterns[snpId]
            
            if( patt ){
                patt.transFac.each{ transFacCounts[it][1]++ }
                if( patt.transFac ){ transFacCounts['allTFs'][1]++ }
                
                patt.motifs.each{
                    motifDists[it.key][1] << it.value
                    motifDists[it.key][3]++ 
                }
                if( patt.motifs ){ 
                    motifDists['allMotifs'][1] += patt.motifs.values()
                    motifDists['allMotifs'][3]++
                }
            }
        }
        
        ////
        // computes fisher exact test using 2 × 2 contingency table:
        //
        // TF present?   | yes | no
        // -----------------------------
        // selected_snps | a   | b
        // -----------------------------
        // all_snps      | c   | d
        //
        
        //create Fisher exact test object
        FisherExact ftest = new FisherExact(patterns.size())
        def results = []
        
        transFacCounts.each{ tf, counts->
            int a = counts[1]
            int b = selectedSnps.size()-a
            int c = counts[0]
            int d = patterns.size()-selectedSnps.size()-c
            
            Double pval = ftest.getCumlativeP(a, b, c, d)
            results << new FisherTestResult( name:tf, grp1:'selected_snps',
                        grp2:'all_snps', a:a, b:b, c:c, d:d, pvalue:pval)
        }
        
        // write output for transFac Fisher tests
        def writer = new PrintWriter(outdir+'transf_'+FeatureSnpsEnrichment.FISHER_TEST_OUT)
        writer.println("name\tgrp1\tgrp2\ta\tb\tc\td\tpvalue") 
        results.each{ writer.println(it.toString()) }
        writer.close()
        
        
        results.clear()
        motifDists.each{ motif, counts->
            int a = counts[3]
            int b = selectedSnps.size()-a
            int c = counts[2]
            int d = patterns.size()-selectedSnps.size()-c
            
            Double pval = ftest.getCumlativeP(a, b, c, d)
            results << new FisherTestResult( name:motif, grp1:'selected_snps',
                        grp2:'all_snps', a:a, b:b, c:c, d:d, pvalue:pval)
        }
        
        // write output for motifs Fisher tests
        writer = new PrintWriter(outdir+'motifs_'+FeatureSnpsEnrichment.FISHER_TEST_OUT)
        writer.println("name\tgrp1\tgrp2\ta\tb\tc\td\tpvalue") 
        results.each{ writer.println(it.toString()) }
        writer.close()
        
        ///////
        // computes one-way Anova for motifs distributions
        writer = new PrintWriter(outdir+ANOVA_TEST_OUT)
        writer.println("name\tgrp1\tmean_grp1\tgrp2\tmean_grp_2\tpvalue")
        
        def avg = { list-> list.isEmpty() ? 'NA' : list.sum()/((double)list.size()) }
        
        motifDists.each{ motif, dist->
            if( [dist[0],dist[1]].any{it.size()<2} ){
                writer.println("${motif}\tall_snps\t${avg(dist[0])}\tselected_sns\t${avg(dist[1])}\tNA")
            }
            else{
                double pval = TestUtils.oneWayAnovaPValue([dist[0] as double[], dist[1] as double[]])
                writer.println("${motif}\tall_snps\t${avg(dist[0])}\tselected_sns\t${avg(dist[1])}\t${pval}")
            }
        }
        
        writer.close()
    }
    
}

