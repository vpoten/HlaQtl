/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import org.biojava.nbio.core.sequence.DNASequence
import org.htmlcleaner.*
import org.msgenetics.hlaqtl.hla.HLAAllele
import org.biojava.nbio.core.sequence.AccessionID
import org.biojava.nbio.core.sequence.template.Sequence
import org.biojava.nbio.core.sequence.compound.DNACompoundSet
import org.biojava.nbio.alignment.Alignments
import org.biojava.bio.program.gff.GFFRecord
import org.biojava.bio.seq.StrandedFeature
import org.ngsutils.Utils
import org.ngsutils.annotation.biojava.SimpleGFFRecordLight
import org.ngsutils.annotation.*

/**
 *
 * @author victor
 */
class HlaSequences {

    static def HLA_ALLELE_URL = 'http://www.ebi.ac.uk/cgi-bin/imgt/hla/get_allele.cgi?'
    static def HLA_EBI_FTP = 'ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/'
    static def HLA_UCSC_DNA = 
        "http://genome.ucsc.edu/cgi-bin/hgc?g=htcGetDna2&c={chr}&l={start}&r={end}" +
        "&db={db}&hgSeq.casing=upper&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA"
    static def HLA_DBFETCH = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=imgthla&id={id}&format=annot&style=raw&Retrieve=Retrieve'
    
    static def HLA_GEN_EXCEPTION = [
        (Main.DRB5):'http://genome.ucsc.edu/cgi-bin/hgc?g=htcDnaNearGene&i=uc003obj.3&db=hg19&c=chr6&l=32485153&r=32498006&o=knownGene&boolshad.hgSeq.promoter=0&boolshad.hgSeq.utrExon5=0&hgSeq.cdsExon=on&boolshad.hgSeq.cdsExon=0&boolshad.hgSeq.utrExon3=0&hgSeq.intron=on&boolshad.hgSeq.intron=0&boolshad.hgSeq.downstream=0&hgSeq.granularity=feature&hgSeq.padding5=0&hgSeq.padding3=0&boolshad.hgSeq.splitCDSUTR=0&hgSeq.casing=exon&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&submit=get+DNA'
        ]
    
    static Map alleleMap = [:]//allele map= key:gene, value:map of {allele_name,HLAAllele}
    static Map positionMap = [:]//pos. map= key:gene, value:map with keys={position,strand}
    static Map<String,Gene> annotMap = [:]//allele map= key:gene, value:annotation object (Gene)
    
    private static intergenic = [:] //map of intergenic DNASequence objects
    
    static String workDir
    static posRegex = /(\w+):(\d+)-(\d+)/
    static int NBASES_LIM = 20000
    
    /**
     * chr6:32403141-32641532
     *
     */
    
    /**
     * loads sequences from external public databases
     */
    static def load(){
        loadPositions()
        Main.hlaGenes.each{ loadGene(it) }
        loadIntergenic()
    }
    
    /**
     * 
     */ 
    static HLAAllele getAllele(String name) {
        def gene = name.substring(0, name.indexOf('*')+1)
        def finalName = getAlleleName(name)
        
        if(!finalName)
            return null
        
        return alleleMap[gene].find{ it.key==finalName }.value
    }
    
    /**
     * gets an allele sequence by name
     * 
     * @return a DNASequence
     */
    static DNASequence getAlleleSeq(String name, boolean reverseIfNeg=true){
        def gene = name.substring(0, name.indexOf('*')+1)
        def allele = getAllele(name)
        
        if(!allele)
            return null
        
        def seq = allele.geneSeq

        if( reverseIfNeg && positionMap[gene]['strand']=='-' )
            return new DNASequence( seq.getReverseComplement().getSequenceAsString() )

        return seq
    }
    
    /**
     * returns the name of the allele that best fit the given name
     */
    static def String getAlleleName(String name, List alleles = null) {
        def gene = name.substring(0, name.indexOf('*')+1)
        def allele = name.substring(name.indexOf('*')+1)
        def toks = allele.split(":")
        
        for(i in toks.length..1){
            allele = gene
            (0..i-1).each{ allele += toks[it]+':' }
            allele = allele.substring(0, allele.length()-1) //remove last ':'
            
            if( alleles==null ) {
                def pair = alleleMap[gene].find{ it.key.startsWith(allele) }
            
                if( pair )
                    return pair.key
            }
            else {
                def alleleObj = alleles.find{ it.name.startsWith(allele) }
                
                if( alleleObj )
                    return alleleObj.name
            }
        }
        
        return null
    }
    
    /**
     * 
     */
    static def writeAnnotation(String name, int offset, String seqId, Writer writer) {
        def allele = getAllele(name)
        
        if(!allele)
            throw new RuntimeException("${name} allele not found")
            
        allele.writeAnnotation(offset, seqId, writer)
    }
    
    
    /**
     * read fasta sequences from input stream
     * 
     * @return a Map<String,DNASequence> or a List<DNASequence>
     */
    public static def readFasta(reader, boolean asList = false){
        def seqs = (asList) ? ([]) : ([:] as TreeMap)
        def seqsIds = [] as Set
        
        def seqId = ''
        def str = null
        
        def createSeq = {
            def dnaSeq = new DNASequence(str.toString(), DNACompoundSet.getDNACompoundSet())
            dnaSeq.accession = new AccessionID(seqId)
            return dnaSeq
        }
        
        reader.eachLine{ line->
            if( !line ) {
                //nothing to do
            }
            else if( line.startsWith('>') ) {
                if( seqId ){
                    seqsIds << seqId
                    if( asList ){ seqs << createSeq() }
                    else{ seqs[seqId] = createSeq() }
                }
                
                seqId = line.substring(1)
                str = new StringBuilder()
            }
            else {
                str << line.trim()
            }
        }
        
        if( !(seqId in seqsIds) ) {//add the last sequence
            if( asList ){ seqs << createSeq() }
            else{ seqs[seqId] = createSeq() }
        }
            
        reader.close()
        return seqs
    }
    
    /**
     *
     */
    private static def readUCSCFasta(URL aUrl, boolean asList = false){
        CleanerProperties props = new CleanerProperties()
        props.setOmitDoctypeDeclaration(true)
        props.setOmitComments(true)
        HtmlCleaner cleaner = new HtmlCleaner(props)
        
        TagNode node = cleaner.clean(aUrl)
        def fasta = node.findElementByName('PRE',true).text
        
        if( !fasta )
            return null
            
        return readFasta( new StringReader(fasta.toString()), asList)
    }
    
    
    /**
     *
     */
//    private static parseHtml(adress){
//        // Clean any messy HTML
//        def cleaner = new HtmlCleaner()
//        def node = cleaner.clean(address.toURL())
//
//        // Convert from HTML to XML
//        def props = cleaner.getProperties()
//        def serializer = new SimpleXmlSerializer(props)
//        def xml = serializer.getXmlAsString(node)
//
//        // Parse the XML into a document we can work with
//        return new XmlParser(false,false).parseText(xml)
//    }
    
    /**
     * 
     * downloads the allele annotation from IMGT and adds it to gene
     */
    private static boolean downloadAnnotation(hlaId, HLAAllele allele, Gene gene, String _geneName) {
        String str = HLA_DBFETCH
        str = str.replace('{id}', hlaId)
        String annotFile = "${hlaId}.unikb"
        
        def reader = createLocalReader(annotFile, str)
        
        // closures for line tokenization
        def extractKey = { line-> line.substring(5,13).trim() }
        def extractDesc = { line-> line.substring(21).trim() }
        
        boolean inExon = false;//flags
        boolean inCDS = false;
        
        def coordRgx = /(\d+)\.\.(\d+)/
        def numberRgx = /\/number="(\d+)"/
        
        SimpleGFFRecordLight current = null;
        
        reader.eachLine{ line->
            if( line.startsWith('FT') ){
                def key = extractKey(line)
                def desc = extractDesc(line)
                
                if( key=='CDS' ){
                    inCDS = true
                }
                else if( key=='exon' ){
                    inExon = true
                    inCDS = false
                    // get coordinates
                    def mat = (desc =~coordRgx)[0]
                    
                    current = new SimpleGFFRecordLight()
                    current.setSeqName( gene.getSeqName() )
                    current.setStart(mat[1] as Integer)
                    current.setEnd(mat[2] as Integer)
                    current.setTranscriptId(allele.name)
                    current.setGeneName(gene.getGeneId())
                    current.setGeneId(gene.getGeneId())
                    ///current.setStrand( 
                    ///    (positionMap[_geneName]['strand']=='+') ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE )
                    current.setStrand(StrandedFeature.POSITIVE) // to simplify, we work in positive strand only
                }
                else if(inExon){
                    inExon = false
                    // get exon number
                    def num = (desc =~numberRgx)[0][1]
                    current.setExonNumber(num)
                    gene.addFeature(current)
                }
                else if(inCDS){
                    if( desc=='/partial' ){
                        allele.isPartial = true
                    }
                }
            }
        }
        
        reader.close()
        
        if( allele.cds.length!=allele.annotationLength ){
            return false
        }
        
        return true
    }
    
    /**
     * 
     */ 
    protected static createLocalReader(name, url) {
        String path = "${workDir}${name}"
        File file = new File(path)
        
        if( !(file.exists()) ) {// download file if not exists
            if( !Utils.download(url, path) ){
                throw new RuntimeException("Cannot download from ${url}")
            }
        }
        
        return Utils.createReader(file)
    }
    
    /**
     * 
     * @return a Map<String,HLAAllele>
     */
    protected static def downloadGene(name){
        ///def root = parseHtml(HLA_ALLELE_URL+hapName)
        ///def seq = root.body.div.find{it.@id=='contents'}.table.tbody.tr.
        ///    td.find{it.@id=='contentsarea'}.table.tbody.tr[25].td.pre.text()
        
        String name2 = ''
        
        if( name in [Main.DRB5, Main.DRB4, Main.DRB3] )//exception DRB345
            name2 = 'DRB345'
        else
            name2 = name.substring(0, name.indexOf('*'))
       
        // load all alleles from *_nuc.fasta file
        def seqs = readFasta( createLocalReader(name2+'_nuc.fasta', HLA_EBI_FTP+name2+'_nuc.fasta') )
        def seqs2 = [:] as TreeMap
        
        Gene gene = new Gene(name.substring(0, name.indexOf('*')), 'HLA')
        annotMap[name] = gene
        
        
        //extract the allele name and add it to the map
        seqs.each{ k,v->
            def toks = k.split("\\s")
            def alleleName = toks[1]
            def hlaId = toks[0].substring(4)
            
            if( alleleName.startsWith(name) ) {
                def allele = new HLAAllele(name:alleleName, hlaId:hlaId, cds:v, gene:gene)
                seqs2[alleleName] = allele
                
                // get annotation for this allele
                if( !downloadAnnotation(hlaId, allele, gene, name) ){
                    // check annotation length != cds length
                    System.err.println("${alleleName} annotation is inconsistent. removed from lists.")
                    seqs2.remove(alleleName)
                }
            }
        }
       
            
        if( name in HLA_GEN_EXCEPTION.keySet() ) {
            //there is no gene sequence in HLA database so download sequence from UCSC
            // the features are in different fasta records (introns in lower-case)
            seqs = readUCSCFasta(HLA_GEN_EXCEPTION[name].toURL(), true)
            
            // build a gene sequence using CDS from HLA database and introns from UCSC
            def intronsIdxs = (1..seqs.size()-1).step(2)
            
            //build geneSeq of parents
            def completeList = seqs2.values().findAll{!it.isPartial}
            completeList*.buildGeneSeqFromIntrons(seqs[intronsIdxs])
            
            seqs2.each{ alleleName, allele->
                if( allele.isPartial ) {
                    //assign parent to partially annotated alleles
                    def match = getAlleleName(alleleName, completeList)
                    allele.parent = seqs2[match]
                    allele.buildGeneSeqFromIntrons(seqs[intronsIdxs])
                }
            }
        }
        else {
            // load alleles with complete sequence from *_gen.fasta file
            if( name in [Main.DQA1, Main.DQB1] )//name exception
                name2 = name.substring(0, name.indexOf('*')-1)
            else
                name2 = name.substring(0, name.indexOf('*'))
            
            
            seqs = readFasta( createLocalReader(name2+'_gen.fasta', HLA_EBI_FTP+name2+'_gen.fasta') )
        
            //extract the allele name and add the gene sequence
            seqs.each{ k,v->
                def toks = k.split("\\s")
                def allele = toks[1]
                seqs2[allele].geneSeq = v
            }
            
            // assign pàrent to each allele
            def completeList = seqs2.values().findAll{it.geneSeq}
            completeList*.rebuildAnnotation()

            seqs2.each{ alleleName, allele->
                if( allele.geneSeq==null ) {
                    def match = getAlleleName(alleleName, completeList)

                    if( !match){
                        System.err.println("No similar match for ${alleleName}, use ${seqs2.firstKey()} as parent")
                        allele.parent = seqs2.firstEntry().value
                    }
                    else{
                        allele.parent = seqs2[match]
                    }
                    
                    allele.buildGeneSeq()
                }
            }//
        }
        
        return seqs2
    }
    
    /**
     * download from UCSC Dna
     * 
     * @return a Map<String,DNASequence>
     */
    protected static def downloadDna(chr, start, end, strand, db = 'hg19'){
        String str = HLA_UCSC_DNA
        str = str.replace('{chr}', chr)
        str = str.replace('{start}', start as String)
        str = str.replace('{end}', end as String)
        str = str.replace('{db}', db)
        
        if( strand=='-') {
            str += '&hgSeq.revComp=on'
        }
        
        def aUrl = str.toURL()
        def seqs = readUCSCFasta(aUrl)
        return seqs
    }
    
    /**
     * loads positions resource file
     */
    protected static loadPositions(){
        def istr = HlaSequences.getClassLoader().getSystemResourceAsStream("hla_pos.txt")
        def reader = new BufferedReader(new InputStreamReader(istr))
        
        reader.splitEachLine("\\s"){ toks->
            positionMap[Main.hlaGenes.find{it.startsWith(toks[0])}] = 
                ['position':toks[1],'strand':toks[2]]
        }
        
        reader.close()
    }
    
    /**
     * 
     */
    protected static loadGene(String gene) {
        alleleMap[gene] = downloadGene(gene)
    }
    
    /**
     * loads intergenic regions
     */
    protected static loadIntergenic(){
        int start = 0
        int end = 0
        def chr = 'chr6'
        
        //load start region
        end = ((positionMap[Main.DRA]['position'] =~ posRegex)[0][2] as Integer) -1
        start = end - NBASES_LIM
        intergenic['_'+Main.DRA] = downloadDna(chr, start, end, '+')
        
        //closure for intergenic region downloading
        def _intergenic = { gene1, gene2->
            start = ((positionMap[gene1]['position'] =~ posRegex)[0][3] as Integer) +1
            end = ((positionMap[gene2]['position'] =~ posRegex)[0][2] as Integer) -1
            intergenic[gene1+'_'+gene2] = downloadDna(chr, start, end, '+')
        }
        
        _intergenic(Main.DRA,Main.DRB5)
        _intergenic(Main.DRB5,Main.DRB1)
        _intergenic(Main.DRB1,Main.DQA1)
        _intergenic(Main.DQA1,Main.DQB1)
        
        //load end region
        start = ((positionMap[Main.DQB1]['position'] =~ posRegex)[0][3] as Integer) +1
        end = start + NBASES_LIM
        intergenic[Main.DQB1+'_'] = downloadDna(chr, start, end, '+')
    }
    
    
    /**
     * gets intergenic sequence between genes
     * 
     * @return a DNASequence
     */
    static def getIntergenic(gene1, gene2){
        if( gene1==Main.DRA && gene2.startsWith('DRB') )
            return intergenic[Main.DRA+'_'+Main.DRB5].firstEntry().value
        else if( gene1.startsWith('DRB') && gene2.startsWith('DRB') )
            return intergenic[Main.DRB5+'_'+Main.DRB1].firstEntry().value
        else if( gene1.startsWith('DRB') && gene2==Main.DQA1 )
            return intergenic[Main.DRB1+'_'+Main.DQA1].firstEntry().value
        else if( !gene1 )
            return intergenic['_'+Main.DRA].firstEntry().value //for testing
        else if( gene2==Main.DQB1 )
            return intergenic[Main.DQA1+'_'+Main.DQB1].firstEntry().value //for testing
            
        return intergenic[gene1+'_'+gene2].firstEntry().value
    }
    
    /**
     * gets the different haplotypes present in a subject list
     * 
     */
    static Set haplotypeSet(subjects){
        def haploset = [] as TreeSet
        
        subjects.each{ subj->
            subj.encodeHaplotypes()
            haploset.addAll( subj.haplotypes )
        }
        
        return haploset
    }
    
    /**
     * gets the different genotypes present in a subject list
     * 
     */
    static Set genotypeSet(subjects){
        return subjects.collect{ it.genotype } as TreeSet
    }
    
    /**
     * @returns a map where key is the genotype and value is a List of subjects 
     * with related genotype
     */
    static Map countGenotypes(subjects) {
        def geno = genotypeSet(subjects)
        def map = [:]
        
        geno.each{ map[it] = [] }
        
        subjects.each{ map[it.genotype] << it }
        
        return map
    }
    
    /**
     * utility method for debugging
     * 
     * @return a map with key=HLAAllele and value=PairwiseSequenceAligner objects
     */
    static def alignAlleles(gene, Sequence target, int total, included) {
        def result = [:]
        def selected = alleleMap[gene].keySet()
        
        if( total!=0 ){
            def alleles = alleleMap[gene].keySet() as List
            selected = [] as Set
            
            while( selected.size()<total ){// select total number of alleles randomly
                selected << alleles[(int)(Math.random()*alleles.size())]
            }
            
            selected += included //add included alleles
        }
        
        selected.each{ name->
            def allObj = alleleMap[gene][name]
            def aligner = allObj.alignment(target, false, false)//global alignment using CDS
            result[allObj] = aligner
        }
        
        Alignments.runPairwiseScorers(result.values() as List)
        
        return result
    }
    
    /**
     * utility method for debugging
     */
    static HLAAllele getBestAllele(mapAligners) {
        def best = mapAligners.keySet().max{
            double sim = mapAligners[it].similarity
            ///println "${it.name} = ${sim}"
            sim
        }
        
        ///println mapAligners[best].pair
        return best
    }
    
    /**
     * utility method for debugging
     */
    static printSequence(gene, name) {
        def allObj = alleleMap[gene][name]
        println allObj.cds.toString()
    }
}

