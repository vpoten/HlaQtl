/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.msgenetics.hlaqtl.Main
import org.ngsutils.Utils
import org.ngsutils.variation.SNPData

/**
 * 
 */ 
class EqtlData {
    String snp
    String isoform
    Double corr
    Double pval
    Double score
}

/**
 * 
 */
class Loci {
    String chr
    int start
    int end
}

/**
 *
 * @author victor
 */
class ImmunobaseAssoc {
	
    //eqtl file fields
    static final int F_SNP = 0
    static final int F_ISO = 1
    static final int F_CORR = 2
    static final int F_PVAL = 3
    static final int F_TYPE = 4
    
    // distances fields
    static final int DIST_EU = 105
    static final int DIST_PW = 106
    static final int DIST_LW = 107
    //
    static public def distEuComp = [ compare:{a,b-> a[DIST_EU]<=>b[DIST_EU]} ] as Comparator
    static public def distPwComp = [ compare:{a,b-> a[DIST_PW]<=>b[DIST_PW]} ] as Comparator
    static public def distLwComp = [ compare:{a,b-> a[DIST_LW]<=>b[DIST_LW]} ] as Comparator
    
    static def distCompMap = [(DIST_EU):distEuComp, (DIST_PW):distPwComp, 
        (DIST_LW):distLwComp]
    
    static def distLblMap = [(DIST_EU):'eucl', (DIST_PW):'pwr', 
        (DIST_LW):'linear']
    
    static final double MISS_PENAL = 0.5*0.5
    static final double MISS_SCORE = 0.1
    static final double MISSING_RATE_MAX = 0.2
    static final double MIN_PVAL = 1e-40
    static final String COMMENT = '#'
    static final String ATT_SNPID = 'marker_mart='
    
    static double EQTL_PVAL_THR = 1e-4
    //LD constants
    double HI_LD_THR = 0.8d
    
    boolean usePval = true //if uses pval or correlation for eqtl score
    def immTabDis = [:] //map with key=disease, value=map of pairs {snp_id,pval}
    def snpsData = [:] as TreeMap //map with key=snp_id, value=SNPData object
    def isoEqtlMap = [:] as TreeMap //map with key=iso_id, value=map of pairs {snp_id,EqtlData}
    def isoGeneMap = [:] as TreeMap //map with key=iso_id, value=gene_name (optional)
    
    //weight function closures
    static def funcLinear = { value-> linearWeight(value, 0.1d, 0.9d) }
    static def funcPower = { value-> powerWeight(value, 2.5d, 0.1d, 0.9d) }
    static def weightFuncMap = [ (DIST_PW):funcPower, (DIST_LW):funcLinear]
    
    List ldSnpsGroups = [[] as TreeSet, [] as TreeSet]
    
    
    
    /**
     * region match for immunobase tab disease data
     */
    Map regionMatch(String locus, String outDir = null) {
        def mat = (locus =~ Utils.locusRegex)[0]
        def chr = mat[1]
        chr = chr.startsWith('chr') ? chr.substring(3) : chr
        int start = mat[2] as Integer
        int end = mat[3] as Integer
        
        //get snps in region
        def regSnps = snpsData.findAll{id, snp-> 
                (snp.chrNum==chr && snp.position>start && snp.position<end)
            }.keySet() as TreeSet
        
        println "\n***** ${regSnps.size()} SNPs in region ${locus}"
        
        // get intersection of eQTL snps and compute the max of eqtl correlation in region
        def isoSnps = [] as TreeSet
        Double eqtlMax = getIntersectSnps(regSnps, isoSnps)
        println "${isoSnps.size()} eQTL SNPs in region"
            
        
        def result = [:] //results: map with key=disease and value=list of isoform scores
        
        immTabDis.keySet().each{ dis->
            // for each disease
            println "\n===== Disease: ${dis}"
            
            // get immunobase snps in region
            def immSnps = immTabDis[dis].findAll{ it.key in regSnps }
            println "${immSnps.size()} immunobase SNPs in region"
            
            normalizeImmunobase(immSnps)
            def weights = immunobaseWeights(immSnps)
            
            // intersection between snps in immunobase and eqtls
            def intersec = immSnps.keySet().findAll{ it in isoSnps }
            println "${intersec.size()} intersection SNPs in region"
            
            result[dis] = calcDistance(immSnps, weights, intersec, eqtlMax)
            
            if( outDir ){
                [DIST_EU, DIST_PW, DIST_LW].each{ distType->
                    def file = "${outDir}${dis}_snps_${locus.replace(':','_').replace('-','_')}_${distLblMap[distType]}.txt"
                    writeRegionDataset(dis, file, intersec, result[dis], distType)
                }
            }
        }
        
        return result
    }
    
    /**
     * normalize immunobase values
     */
    protected void normalizeImmunobase(Map immSnps) {
        Double invmax = immSnps.max{it.value}?.value
        invmax = (invmax!=null) ? 1.0/invmax : null
        immSnps.keySet().each{ immSnps[it] = immSnps[it]*invmax }
    }
    
    /**
     * immunobase Weights internal
     */ 
    protected _immunobaseWeights(Map normImmSnps, function) {
        def weights = [:] as TreeMap
        double sum = 0.0
        
        normImmSnps.each{ snpId, value->
            double weight = function(value)
            sum += weight
            weights[snpId] = weight
        }
        
        sum = 1.0/sum
        weights.keySet().each{ weights[it] = weights[it]*sum }
        
        return weights
    }
    
    /**
     * 
     */
    protected immunobaseWeights(Map normImmSnps) {
        def weights = [:]
        weightFuncMap.each{ dist, function->
            weights[dist] = _immunobaseWeights(normImmSnps, function)
        }
        return weights
    }
    
    /**
     * function
     */
    static protected double linearWeight(double value, double min, double max){
        return (max-min)*value + min
    }
     
    /**
     * function
     */
    static protected double powerWeight(double value, double exp, double min, double max){
        return (max-min)*Math.pow(value, exp) + min
    }
    
    /**
     * intersects isoEqtlMap snps vs snpSet
     * 
     * @return the max score for eqtls with snps in snpSet
     */ 
    protected double getIntersectSnps(Set snpSet, Set intersec) {
        Double eqtlMax = Double.MIN_VALUE
        isoEqtlMap.each{ isoId, isoEqtls->
            isoEqtls.each{ 
                if( it.key in snpSet ){ 
                    intersec << it.key
                    if( it.value.score>eqtlMax ){ eqtlMax = it.value.score }
                } 
            }
        }
        
        return eqtlMax
    }
    
    /**
     * 
     */ 
    protected List calcDistance(Map immSnps, Map weights, Set intersec, Double eqtlMax) {
        def result = []
        
        isoEqtlMap.each{ isoId, isoEqtls->
            // calc. distance between isoform eqtls and immunobase snps
            if( missingRate(isoEqtls, intersec)<MISSING_RATE_MAX ){ // exclude isoforms with high rate of missing snps
                double dist = calcEuclDistance(immSnps, isoEqtls, null/*eqtlMax*/, intersec)
                def isoResult = [(F_ISO):isoId, (DIST_EU):dist]

                weights.each{ distType, distWeights->
                    dist = weightedDistance(immSnps, isoEqtls, null/*eqtlMax*/, intersec, distWeights)
                    isoResult[distType] = dist
                }

                result << isoResult
            }
        }
        
        return result
    }
    
    /**
     * 
     */ 
    protected loadSubjects() {
        if( !Main.allSubjects ){
            Main.loadSubjects(['--group=CEU,GBR,FIN,TSI'])
        }
    }
    
    /**
     * gff3 vs isoforms eqtl match
     */
    Map gff3Match(fileGff, disease) {
        def immSnps = parseGff3(fileGff)// get immunobase snps:score
        
        println "\n===== ${disease}"
        println " ${immSnps.size()} SNPs in ${fileGff}"
        
        
        normalizeImmunobase(immSnps)
        def weights = immunobaseWeights(immSnps)
        
            
        // get eQTL snps in gff3 file and compute the max of eqtl correlation in region
        def intersec = [] as TreeSet
        Double eqtlMax = getIntersectSnps(immSnps.keySet(), intersec)
        println "${intersec.size()} eQTL SNPs in region"
        
        def result = [:] //results: map with key=disease and value=list of isoform scores
        result[disease] = calcDistance(immSnps, weights, intersec, eqtlMax)
            
        return result
    }
    
    /**
     * get snps with LD rSq > HI_LD_THR
     */
    protected getHiLDSnps(bestSnp, resultLD) {
        def hiLdSnps = null
        double ldThr = HI_LD_THR
        while( !hiLdSnps ) { 
            hiLdSnps = resultLD[bestSnp].findAll{it.value>ldThr}
            if(hiLdSnps){ 
                println "Used ${ldThr} LD rSquare threshold" 
                println "${hiLdSnps.size()} snps remaining after LD filtering" 
            }
            ldThr -= 0.1
        }
        return hiLdSnps.keySet()+bestSnp
    }
    
    /**
     * 
     */ 
    Map gff3MatchLD(fileGff, String disease, String locus, String outDir) {
        def immSnps = parseGff3(fileGff)// get immunobase snps:score
        
        println "\n===== ${disease}"
        println " ${immSnps.size()} SNPs in ${fileGff}"
        
        // get eQTL snps in gff3 file
        def intersec = [] as TreeSet
        Double eqtlMax = getIntersectSnps(immSnps.keySet(), intersec)
        println "${intersec.size()} eQTL SNPs in region"
        
        // get best immunobase snp
        def bestSnp = immSnps.findAll{it.key in intersec}.max{it.value}.key
        println "Best snp ${bestSnp} with score ${immSnps[bestSnp]}"
        
        normalizeImmunobase(immSnps)
        def weights = immunobaseWeights(immSnps)

        // get subjects
        loadSubjects()
        
        // compute LD
        def resultLD = LDCalculator.perform(locus, [bestSnp], intersec, Main.allSubjects.keySet(), outDir)
        
        // get snps with LD rSq > HI_LD_THR
        def hiLdSnpsSet = getHiLDSnps(bestSnp, resultLD)
        
        // use high LD snps to compute distances
        def result = [:] //results: map with key=disease and value=list of isoform scores
        result[disease] = calcDistance(immSnps, weights, hiLdSnpsSet, eqtlMax)
            
        return result
    }
    
    
    /**
     * add snps in gff3 file to LD groups (calculates LD using 1000genomes data)
     */
    def addGff3ToLDGroups(fileGff, String locus, String outDir) {
        // get immunobase snps:score and populate snpData
        def immSnps = parseGff3(fileGff, [:] as TreeMap, true)
        
        println "\n${immSnps.size()} SNPs in ${fileGff}"
        
        // get eQTL snps in gff3 file
        def intersec = [] as TreeSet
        Double eqtlMax = getIntersectSnps(immSnps.keySet(), intersec)
        println "${intersec.size()} eQTL SNPs in region"
        
        // get best immunobase snp
        def bestSnp = immSnps.findAll{it.key in intersec}.max{it.value}.key
        println "Best snp ${bestSnp} with score ${immSnps[bestSnp]}"
        
        // get subjects
        loadSubjects()
        
        // compute LD
        def resultLD = LDCalculator.perform(locus, [bestSnp], intersec, Main.allSubjects.keySet(), outDir)
        
        // get snps with LD rSq > HI_LD_THR
        def hiLdSnpsSet = getHiLDSnps(bestSnp, resultLD)
        
        ldSnpsGroups[0] += intersec.findAll{it in hiLdSnpsSet}
        ldSnpsGroups[1] += intersec.findAll{!(it in hiLdSnpsSet)}
    }
    
    /**
     * 
     */ 
    def writeLdSnpsGroups(file) {
        def writer = new File(file).newPrintWriter() 
        
        //print header
        def header = "id\tchr\tposition\talleles\tminor\tmajor\tmaf\tattributes"
        writer.println(header)
        
        ldSnpsGroups.sum().each{ snp->
            def data = snpsData[snp]
            def group = (snp in ldSnpsGroups[0]) ? 'high' : 'low'
            writer.println(
"${snp}\t${data.chr}\t${data.position}\t${data.alleles}\t${data.minor}\t${data.major}\t${data.maf}\t${group}"
            )
        }
        
        writer.close()
    }
    
    /**
     *
     */
    def printResult(resultMap, writer, String locus, int num, int distType = DIST_EU){
        resultMap.each{ dis, list->
            // sort isoforms by score/distance (asc)
            Collections.sort(list, distCompMap[distType])
        
            (0..num-1).each{
                if( it<list.size() ){
                    def isoform = list[it][F_ISO]
                    def gene = isoGeneMap[isoform]
                    writer.println("${dis}\t${locus}\t${isoform}\t${list[it][distType]}\t${gene?:'NA'}")
                }
            }
        }
    }
    
    /**
     *
     */
    def printResultHeader(writer){
        String header = "disease\tregion\tisoform\tdistance\tgene"
        writer.println(header)
    }
    
    /**
     *
     */
    protected def writeRegionDataset(disease, file, snpSet, scores, int distType = DIST_EU) {
        // sort isoforms by score/distance (asc)
        Collections.sort(scores, distCompMap[distType])
            
        //sort snps by position
        def snpList = snpSet.collect{ snpsData[it] }
        Collections.sort(snpList, SNPData.comparator)
        
        //select isoforms with enough markers
        double scoreLimit = Math.sqrt(MISS_PENAL*0.3d*snpSet.size())
        def isoforms = scores.findAll{ it[distType]<scoreLimit }.collect{it[F_ISO]}
        
        def formatScore = { val-> (val==null) ? 'NA' : val }
        
        //print header
        def header = "snp\tlocus\tmaf\t${disease}"+((isoforms) ? isoforms.sum{"\t${it}"} : '')
        
        def writer = new File(file).newPrintWriter() 
        writer.println(header)
        
        //print snp info and scores for chip and eqtls
        snpList.each{ snp->
            writer.print("${snp.id}\t${snp.locus}\t${snp.maf}\t${formatScore(immTabDis[disease][snp.id])}")
            writer.println( (isoforms) ? isoforms.sum{"\t${formatScore(isoEqtlMap[it][snp.id]?.score)}"} : '' )
        }
        
        writer.close()
    }
    
    /**
     *
     */
    protected double missingRate(isoEqtls, snpsSet) {
        double total = snpsSet.sum{ isoEqtls.containsKey(it) ? 0.0d : 1.0d }
        return total/(double)snpsSet.size()
    }
    
    /**
     *
     */
    protected double isoEqtlsMaxScore(isoEqtls, snpsSet) {
        Double eqtlMax = Double.MIN_VALUE
        
        isoEqtls.each{ 
            if( it.key in snpsSet ){
                if( it.value.score>eqtlMax ){ eqtlMax = it.value.score }
            } 
        }
        
        return eqtlMax
    }
    
    /**
     *
     */
    protected double calcEuclDistance(immSnps, isoEqtls, Double isoMax, snpsSet){
        if( isoMax==null ) {
            isoMax = isoEqtlsMaxScore(isoEqtls, snpsSet)
        }
        
        double dist = 0.0
        isoMax = 1.0/isoMax

        snpsSet.each{
            Double x1 = immSnps[it]
            Double x2 = isoEqtls[it]?.score
            
            x2 = (x2==null) ? MISS_SCORE : x2*isoMax
            
            double val = x1-x2
            dist += val*val
        }
        
        return Math.sqrt(dist)/((double)snpsSet.size())//normalize distance
    }
    
    /**
     *
     */
    protected double weightedDistance(immSnps, isoEqtls, Double isoMax, snpsSet, weights){
        if( isoMax==null ) {
            isoMax = isoEqtlsMaxScore(isoEqtls, snpsSet)
        }
        
        double dist = 0.0
        isoMax = 1.0/isoMax

        snpsSet.each{
            Double x1 = immSnps[it]
            Double x2 = isoEqtls[it]?.score
            double w = weights[it]
            
            x2 = (x2==null) ? MISS_SCORE : x2*isoMax
            
            dist += w*Math.abs(x1-x2)
        }
        
        return dist/((double)snpsSet.size())//normalize distance
    }
    
    
    /**
     *
     */
    def parseImmunobaseTab(file) {
        def reader = Utils.createReader(new File(file))
        def header = reader.readLine()
        
        def fields = header.split("\t")
        int ndiseases = fields.length-2
        
        (0..ndiseases-1).each{ immTabDis[fields[2+it]] = [:] as TreeMap }
        
        reader.splitEachLine("\t"){ toks->
            int id = toks[1] as Integer
            String strId = (id>0) ? "rs${id}" : id
            
            (0..ndiseases-1).each{
                def valStr = toks[2+it]
                
                if( valStr!='NA' ) {
                    valStr = valStr.startsWith('>') ? valStr.substring(1) : valStr
                    double val = valStr as Double
                    val = (val<MIN_PVAL) ? MIN_PVAL : val // saturate pvalues
                    val = -Math.log10(val) //transform val to log scale
                    immTabDis[fields[2+it]][strId] = val
                }
            }
        }
        
        reader.close()
    }
    
    /**
     * 
     *  @return a map with key=snpId and value=score
     */ 
    public Map parseGff3(file, map = [:] as TreeMap, boolean fillSnpData = false) {
        def reader = Utils.createReader(new File(file))
        def header = reader.readLine()
        
        reader.splitEachLine("\t"){ toks->
            if( !toks[0].startsWith(COMMENT) ){
                def attributes = toks[8].split(';') as List
                def id = attributes.find{it.startsWith(ATT_SNPID)}?.substring(ATT_SNPID.length())
                String strId = (id) ? "rs${id}" : null

                if( strId ) {
                    map[strId] = (toks[5] as Double)//score is already in log scale
                    
                    if( fillSnpData && !snpsData.containsKey(strId) ){
                        //add snp to snpData
                        snpsData[strId] = new SNPData( id:strId, chr:toks[0], 
                            position:toks[3] as Integer,
                            minor:'?', major:'?', maf:0.0d )
                    }
                }
            }
        }
        
        reader.close()
        
        return map
    }
    
    /**
     *
     */
    def parseImmunobaseGff3(file, disease) {
        def reader = Utils.createReader(new File(file))
        def header = reader.readLine()
        
        def mapGff = parseGff3(file)
        
        def map = immTabDis[disease]
        
        if( map==null ){
            map = [:] as TreeMap
            immTabDis[disease] = map
        }
        
        map.putAll(mapGff)
    }
    
    /**
     * add snps to internal SNPData map
     */
    def parseSnpsLog(snpsFile){
       // parse snps log file
        def reader = Utils.createReader(new File(snpsFile))
        reader.readLine() //skip header
        
        reader.splitEachLine("\t"){ toks->
            def snp = new SNPData( 
                id:toks[0], chr:toks[1], position:toks[2] as Integer,
                minor:toks[4], major:toks[5], maf:toks[6] as Double )
            
            //add snp to snpData
            snpsData[toks[0]] = snp
        }
        
        reader.close() 
    }
    
    /**
     * 
     * add eqtls to internal isoform EqtlData map
     * @param file : eqtl file
     * @param pPval : eqtl pvalue threshold
     * @param loci : list of desired loci (eqtl snps inside them) or null to disable filtering
     */
    def parseEqtlFile(file, double pPval=0.05, loci=null) {
        def listLoci = null
        def missSnps = null
        
        if( loci ) {
            listLoci = []
            missSnps = [] as TreeSet
            
            loci.each{ locus->
                def mat = (locus =~ Utils.locusRegex)[0]
                def chr = mat[1]
                chr = chr.startsWith('chr') ? chr.substring(3) : chr
                int start = mat[2] as Integer
                int end = mat[3] as Integer
                
                listLoci << new Loci(chr:chr, start:start, end:end)
            }
        }//
        
        def inLoci = { snpId->
            if( !listLoci ){ return true }
            def snp = snpsData[snpId]
            if( !snp ){ missSnps << snpId; return false }
            listLoci.any{ (snp.chrNum==it.chr && snp.position>it.start && snp.position<it.end) }
        }//
        
        parseEqtlFileInternal(file, pPval, inLoci, isoEqtlMap)
        
        if( missSnps ) {
            println "parseEqtlFile: Warning ${missSnps.size()} missing snps."
        }
    }
    
    /**
     * 
     * add eqtls to internal isoform EqtlData map
     * @param file : eqtl file
     * @param pPval : eqtl pvalue threshold
     * @param snpSet : list of allowed eqtl snps
     */
    def parseEqtlFileBySnps(file, double pPval=0.05, snpsSet) {
        def filter = { snpId-> return (snpId in snpsSet) }
        parseEqtlFileInternal(file, pPval, filter, isoEqtlMap)
    }
    
    /**
     *
     */
    protected def parseEqtlFileInternal(file, double pPval, filter, eqtlMap) {
        def reader = Utils.createReader(new File(file))
        reader.readLine() //skip header
        
        reader.splitEachLine("\t"){ toks->
            //fields: snp isoform correlation pvalue type
            double pval = toks[F_PVAL] as Double
            def snpId = toks[F_SNP]

            if( pval<pPval && toks[F_TYPE]!='trans' && filter(snpId) ){
                
                def isoId = toks[F_ISO]
                
                //store isoform eqtl
                if( eqtlMap[isoId]==null )
                    eqtlMap[isoId] = [:] as TreeMap
                    
                double corr = Math.abs(toks[F_CORR] as Double) 
                double score = (usePval) ? -Math.log10(pval) : corr

                eqtlMap[isoId][snpId] = new EqtlData( isoform:isoId, 
                    snp:snpId, corr:corr, pval:pval, score:score )
            }
        }
        
        reader.close()
        
        // remove isoforms with non significant eqtls
        def toRemove = eqtlMap.keySet().findAll{ isoId-> !eqtlMap[isoId].any{it.value.pval<EQTL_PVAL_THR} }
        toRemove.each{ eqtlMap.remove(it) }
        
    }
    
}

