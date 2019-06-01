package org.msgenetics.hlaqtl

import au.com.bytecode.opencsv.CSVReader
import org.biojava.nbio.core.sequence.io.FastaWriterHelper
import org.biojava.nbio.core.sequence.DNASequence
import org.biojava.nbio.core.sequence.AccessionID
import org.ngsutils.Utils
import org.ngsutils.AnnotationDB
import org.awsjoblauncher.LocalJobSlave
import org.awsjoblauncher.SqsJobMaster
import org.awsjoblauncher.SqsJobSlave
import org.awsjoblauncher.SqsAuditReceiver
import org.awsjoblauncher.SqsAuditSender
import org.awsjoblauncher.postprocess.PostBase
import org.awsjoblauncher.postprocess.JobResponse
import org.awsjoblauncher.postprocess.CommConstants
import org.awsjoblauncher.storage.*
import org.awsjoblauncher.QueueUtils
import org.msgenetics.hlaqtl.eqtl.*
import org.ngsutils.eqtl.CorrelationCalc


/**
 * HlaQtl Main class
 *
 */
public class Main 
{
    
    //gene keys
    static def DRA = 'DRA*'
    static def DRB1 = 'DRB1*'
    static def DRB5 = 'DRB5*'
    static def DRB4 = 'DRB4*'
    static def DRB3 = 'DRB3*'
    static def DQA1 = 'DQA1*'
    static def DQB1 = 'DQB1*'
    
    static def hlaGenes = [DRA, DRB1, DRB5, DRB4, DRB3, DQA1, DQB1]
    
    static Map subjects = [:] as TreeMap //map of Subject objects
    static Map allSubjects = [:] as TreeMap //map of Subject objects (every geuvadis subject)
    
    static def config //config properties object
    
    //commandline options
    public static final String OPT_OUTPUT_DIR = "--output="
    public static final String OPT_PROPS = "--props="
    public static final String OPT_SUBJECTS = "--subjs="
    public static final String OPT_START_PHASE = "--start="
    public static final String OPT_LOCUS = "--locus="
    public static final String OPT_VCF = "--vcf="
    public static final String OPT_GROUP = "--group="
    public static final String OPT_CUFFPRE = "--cuffPre="
    public static final String OPT_S3OUT = "--s3Out="
    public static final String OPT_EXPRDATA = "--exprData="
    public static final String OPT_TRANSEQTLSNP = "--transEqtlSnps="
    public static final String OPT_PAIRS = "--pairs="
    public static final String OPT_INPUT_DIR = "--input="
    public static final String OPT_SNPS = "--snps="
    public static final String OPT_TAXID = "--taxId="
    public static final String OPT_VALUE = "--value="
    public static final String OPT_BLAST = "--blast="
    public static final String OPT_EQTL_DIR = '--eqtlDir='
    public static final String OPT_ISOFORMS = '--isoforms='
    
    public static final String COMM_QC = 'qc'
    public static final String COMM_FILTER = 'filter'
    public static final String COMM_HLAMAP = 'hlaMap'
    public static final String COMM_CBCBMAP = 'cbcbMap'
    public static final String COMM_EQTL = 'eqtl'
    public static final String COMM_ASSOCORR = 'assocCorr'
    public static final String COMM_DEBUGEQTL = 'debugEqtl'
    public static final String COMM_DOWN_EXPR = 'downExprData'
    public static final String COMM_DOWN_JOBDATA = 'downJobData'
    public static final String COMM_DOWN_EQTLRES = 'downEqtlResult'
    public static final String COMM_PATH_PATTERNS = 'pathwayPatterns'
    public static final String COMM_SNP_GENE = 'snpGeneRel'
    public static final String COMM_EQTLTRANS_FILTER = 'eqtlTransFilter'
    public static final String COMM_DOWN_SEQS = 'downSeqs'
    public static final String COMM_SEARCH_CT = 'searchCisTrans'
    public static final String COMM_ASSOC_EQTL_SNPS = 'assocEqtlSnps'
    public static final String COMM_GO_ENRICH = 'goEnrichment'
    public static final String COMM_GROUP_LD = 'groupLD'
    public static final String COMM_SIMUL_ERR = 'simulErrCorr'
    public static final String COMM_FEATSNP_ENRICH = 'featSnpsEnrichment'
    public static final String COMM_ADDFEAT_SNP = 'addFeatToSnps'
    public static final String COMM_EQTL_TRANS_LD = 'eqtlTransLD'
    public static final String COMM_CNTRL_GROUP_LD = 'controlGroupLD'
    public static final String COMM_GET_GENO = 'getGenotypes'
    public static final String COMM_EQTL_SIMPLE = 'eqtlSimple'
    public static final String COMM_DIFF_EXPR_PIPE = 'diffExprPipe'
    public static final String COMM_PHLAT = 'phlat'
    public static final String COMM_GTEX_SEARCH = 'gtexSearch'
    
    public static def COMMANDS = [COMM_QC, COMM_FILTER, COMM_HLAMAP, COMM_CBCBMAP, 
        COMM_EQTL, COMM_ASSOCORR, COMM_DEBUGEQTL, COMM_DOWN_EXPR, COMM_DOWN_JOBDATA,
        COMM_DOWN_EQTLRES, COMM_PATH_PATTERNS, COMM_SNP_GENE, COMM_EQTLTRANS_FILTER,
        COMM_DOWN_SEQS, COMM_SEARCH_CT, COMM_ASSOC_EQTL_SNPS, COMM_GO_ENRICH,
        COMM_GROUP_LD, COMM_SIMUL_ERR, COMM_FEATSNP_ENRICH, COMM_ADDFEAT_SNP,
        COMM_EQTL_TRANS_LD, COMM_CNTRL_GROUP_LD, COMM_GET_GENO, COMM_EQTL_SIMPLE,
        COMM_DIFF_EXPR_PIPE, COMM_PHLAT, COMM_GTEX_SEARCH]
    
    /** List of commands that skip the pre-load of properties and subjects */
    public static def SKIP_STARTS_LOADS = [COMM_GTEX_SEARCH]
    
    //usages
    public static final String USA_QC = "${COMM_QC} ${OPT_PROPS}<ngsengine props file> [${OPT_SUBJECTS}=]"
    
    public static final String USA_FILTER = "${COMM_FILTER} ${OPT_PROPS}<ngsengine props file> [${OPT_SUBJECTS}=]"
    
    public static final String USA_HLAMAP = "${COMM_HLAMAP} ${OPT_PROPS}<ngsengine props file> " + 
        "${OPT_OUTPUT_DIR}<output_dir> ${OPT_INPUT_DIR}<fastx toolkit base> ${OPT_SUBJECTS}<csv file>"
    
    public static final String USA_CBCBMAP = "${COMM_CBCBMAP} ${OPT_PROPS}<ngsengine props file> " + 
        "[${OPT_SUBJECTS}=] [${OPT_START_PHASE}=]"
    
    public static final String USA_EQTL = "${COMM_EQTL} ${OPT_PROPS}<ngsengine props file> [${OPT_SUBJECTS}=] "+
        "${OPT_OUTPUT_DIR}<output_dir> ${OPT_LOCUS}<chr>:<start>:<end> ${OPT_VCF}<vcf_file> "+
        "${OPT_GROUP}<group> ${OPT_CUFFPRE}=<cufflinks out prefix> [${OPT_EXPRDATA}<dir>] "+
        "[${OPT_TRANSEQTLSNP}<snps file>] [${OPT_ISOFORMS}<isoforms filter file>]"
    
    public static final String USA_ASSOCORR = "${COMM_ASSOCORR} ${OPT_OUTPUT_DIR}<output_dir>"
    
    public static final String USA_DOWN_EXPR = "${COMM_DOWN_EXPR} ${OPT_PROPS}<ngsengine props file> " + 
        "${OPT_OUTPUT_DIR}<output_dir> ${OPT_CUFFPRE}=<cufflinks out prefix> ${OPT_GROUP}<group>"
    
    public static final String USA_DOWN_JOBDATA = "${COMM_DOWN_JOBDATA} ${OPT_PROPS}<ngsengine props file> " + 
        "${OPT_OUTPUT_DIR}<output_dir> ${OPT_GROUP}<group>"
    
    public static final String USA_DOWN_EQTLRES = "${COMM_DOWN_EQTLRES} ${OPT_PROPS}<ngsengine props file> " +
        "${OPT_LOCUS}<chr>:<start>:<end> {OPT_GROUP}<group> ${OPT_CUFFPRE}=<cufflinks out prefix>"
    
    public static final String USA_PATH_PATTERNS = "${COMM_PATH_PATTERNS} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_SUBJECTS}<risk snps arff file> ${OPT_PAIRS}<pairs text file SNP-Gene> ${OPT_TAXID}<ncbi taxonomy id>"
    
    public static final String USA_SNP_GENE = "${COMM_SNP_GENE} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_INPUT_DIR}<eqtl results dir> ${OPT_SNPS}<snps file> ${OPT_TAXID}<ncbi taxonomy id> "+
        "[${OPT_PAIRS}<pairs text file SNP-Gene>]"
    
    public static final String USA_EQTL_STATS = "${COMM_EQTLTRANS_FILTER} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_INPUT_DIR}<eqtl results dir> ${OPT_VALUE}<filter percentile [0,100]>"+
        "${OPT_BLAST}=<blast+ bin path> ${OPT_SNPS}=<cis eqtls file> ${OPT_EQTL_DIR}<eqtl-cis results dir>"
    
    public static final String USA_DOWN_SEQS = "${COMM_DOWN_SEQS} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_VALUE}=<isoforms file>"
    
    public static final String USA_SEARCH_CT = "${COMM_SEARCH_CT} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_INPUT_DIR}<eqtl trans results dir> ${OPT_EQTL_DIR}<eqtl-cis results dir>"
    
    public static final String USA_ASSOC_EQTL_SNPS = "${COMM_ASSOC_EQTL_SNPS} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_INPUT_DIR}<eqtl trans results dir> ${OPT_EQTL_DIR}<eqtl-cis results dir> ${OPT_SNPS}<cis eqtls file> "+
        "${OPT_VALUE}<featured snps>"
    
    public static final String USA_GO_ENRICH = "${COMM_GO_ENRICH} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_PAIRS}<pairs text file feature-set> ${OPT_TAXID}<ncbi taxonomy id>"
    
    public static final String USA_GROUP_LD = "${COMM_GROUP_LD} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_EQTL_DIR}<eqtl-cis results dir> ${OPT_SNPS}<groups file>"
    
    public static final String USA_SIMUL_ERR = "${COMM_SIMUL_ERR} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_EQTL_DIR}<eqtl-cis results dir> ${OPT_ISOFORMS}<filtered isoforms>"
    
    public static final String USA_FEATSNP_ENRICH = "${COMM_FEATSNP_ENRICH} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_SNPS}<snps log group file>"
    
    public static final String USA_ADDFEAT_SNP = "${COMM_ADDFEAT_SNP} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_SNPS}<snps log file>"
    
    public static final String USA_EQTL_TRANS_LD = "${COMM_EQTL_TRANS_LD} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_INPUT_DIR}<eqtl trans results dir> ${OPT_EQTL_DIR}<eqtl-cis results dir> "+
        "${OPT_ISOFORMS}<isoforms filtered file>"
    
    public static final String USA_CNTRL_GROUP_LD = "${COMM_CNTRL_GROUP_LD} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_EQTL_DIR}<eqtl-cis results dir> ${OPT_SNPS}<groups file>"
    
    public static final String USA_GET_GENO = "${COMM_GET_GENO} ${OPT_OUTPUT_DIR}<output_dir> " +
        "${OPT_EQTL_DIR}<eqtl-cis results dir> ${OPT_SNPS}<groups file>"
    
    public static final String USA_EQTL_SIMPLE = "${COMM_EQTL_SIMPLE} ${OPT_PROPS}<ngsengine props file> "+
        "${OPT_OUTPUT_DIR}<output_dir> ${OPT_LOCUS}<chr>:<start>:<end> [${OPT_GROUP}<group> ${OPT_VALUE}hapmap] "+
        "${OPT_EXPRDATA}<file>"
    
    public static final String USA_DIFF_EXPR_PIPE = "${COMM_DIFF_EXPR_PIPE} ${OPT_PROPS}<ngsengine props file> "+\
        "${OPT_SUBJECTS}=<experiments desc. file> ${OPT_START_PHASE}<phase> ${OPT_OUTPUT_DIR}<output_dir>"
    
    public static final String USA_PHLAT = "${COMM_PHLAT} ${OPT_PROPS}<ngsengine props file> "+\
        "${OPT_SUBJECTS}=<subjects file> ${OPT_OUTPUT_DIR}<output_dir>"
    
    //sqs master
    static SqsJobMaster sqsMaster
    static SqsJobSlave sqsSlave
    static SqsAuditSender sqsSender
    static SqsAuditReceiver sqsReceiver
    
    //finished job set
    static def finishedJobs = Collections.synchronizedSet(new HashSet())
    
    //constants for file results
    public static final String CORREL_FILE = 'correlation.txt'
    public static final String MERG_GFF_FILE = 'iso_merged.gtf'
    public static final String MERG_TRACK_FILE = 'iso_merged.tracking'
    public static final String LINE_PRE = '===== '
    public static final int MIN_OBS = 1
    public static final int CORR_FILE_CORR = 2
    public static final int CORR_FILE_PVAL = 3
    public static final int CORR_FILE_CPVAL = 5
    public static final int CORR_FILE_NOBS = 7
    public static final String CORR_DEBUG_FILE = 'debug_corrData.txt'
    public static final String GEUV_FILE = 'E-GEUV-1.sdrf.txt'
    
    public static final String _1000G_DATA = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/'
    
    /**
     * 
     */
    public static void main( String[] args )
    {
        if( args.length==0 || !(args[0] in COMMANDS) ) {
            // first check args
            if( args.length>0 )
                System.err.println( args[0]+" operation, not valid.");

            System.err.println( "Valid operations : "+COMMANDS )
            System.exit(1)
        }
        
        println "Start time: ${new Date()}\n"
        
        if(!(args[0] in SKIP_STARTS_LOADS)) {
            loadConfig(args as List)
            loadSubjects(args as List)
        }
        
        if( args[0] == COMM_QC ){
            launchQC(args as List)
        }
        else if( args[0] == COMM_FILTER ){
            launchFiltering(args as List)
        }
        else if( args[0] == COMM_HLAMAP ){
            hlaMap(args as List)
        }
        else if( args[0] == COMM_CBCBMAP ){
            launchCBCB(args as List)
        }
        else if( args[0] == COMM_EQTL ){
            launchEqtl(args as List)
        }
        else if( args[0] == COMM_ASSOCORR ){
            AssocCorr.perform(args as List)
        }
        else if( args[0] == COMM_DEBUGEQTL ){
            DebugEqtl.perform(args as List, allSubjects)
        }
        else if( args[0] == COMM_DOWN_EXPR ){
            downExpresData(args as List)
        }
        else if( args[0] == COMM_DOWN_JOBDATA ){
            DownJobData.downloadJobData(args as List, allSubjects)
        }
        else if( args[0] == COMM_DOWN_EQTLRES ){
            DownJobData.downloadEqtlResult(args as List)
        }
        else if( args[0] == COMM_PATH_PATTERNS ){
            PathwayPatterns.perform(args as List)
        }
        else if( args[0] == COMM_SNP_GENE ){
            SnpGeneRelation.perform(args as List)
        }
        else if( args[0] == COMM_EQTLTRANS_FILTER ){
            EqtlStatistics.perform(args as List)
        }
        else if( args[0] == COMM_DOWN_SEQS ){
            DownSequences.perform(args as List)
        }
        else if( args[0] == COMM_SEARCH_CT ){
            EqtlStatistics.performSearchCisTrans(args as List)
        }
        else if( args[0] == COMM_ASSOC_EQTL_SNPS ){
            EqtlStatistics.performAssocEqtlSnps(args as List)
        }
        else if( args[0] == COMM_GO_ENRICH ){
            GOEnrichment.perform(args as List)
        }
        else if( args[0] == COMM_GROUP_LD ){
            EqtlSearcher.performGroupLD(args as List)
        }
        else if( args[0] == COMM_SIMUL_ERR ){
            ErrorSimul.perform(args as List, config.org.webngs.tool.threadCount as Integer)
        }
        else if( args[0] == COMM_FEATSNP_ENRICH ){
           FeatureSnpsEnrichment.perform(args as List)
        }
        else if( args[0] == COMM_ADDFEAT_SNP ){
            FeatureSnpsEnrichment.addFeaturesToSnp(args as List)
        }
        else if( args[0] == COMM_EQTL_TRANS_LD ){
            EqtlSearcher.performActualEqtlTransLD(args as List)
        }
        else if( args[0] == COMM_CNTRL_GROUP_LD ){
            EqtlSearcher.performControlGroupLD(args as List)
        }
        else if( args[0] == COMM_GET_GENO ){
            EqtlSearcher.generateGenotypes(args as List)
        }
        else if( args[0] == COMM_EQTL_SIMPLE ){
            EqtlSimpleCalc.performCalc(args as List)
        }
        else if( args[0] == COMM_DIFF_EXPR_PIPE ){
            DiffExprPipeline.perform(args as List)
        }
        else if( args[0] == COMM_PHLAT ){
            phlatTyping(args as List)
        }
        else if( args[0] == COMM_GTEX_SEARCH ){
            GTExSearcher.main(args)
        }
        
        println "End time: ${new Date()}\n"
    }
    
    /**
     * checks that the given dirs exists
     * 
     * @param list : a list of string options
     * @return a list of dir names or null if error
     */
    public static checkAndCleanDirOpts(list){
        def dirs = list.collect{ //extract dir names from options
                if(it){
                    String dir = it.substring( it.indexOf('=')+1 )
                    //add '/' to the dir path
                    dir.endsWith(File.separator) ? dir : dir+File.separator
                }
                else{ '' }
            }
        
        for( dir in dirs ){ //check if the directories exists
            if( !dir || !(new File(dir).exists()) ){
                System.err.println( "directory option missing or not exists.");
                return null
            }
        }
        
        return dirs
    }
    
    /**
     *
     */
    public static def readOption(args, opt){
        def str = args.find{ it.startsWith(opt) }
        
        if( !str  )
            return null
            
        return str.substring(str.indexOf('=')+1)
    }
    
    /**
     *
     */
    public static Properties getProperties(){
        return config.toProperties()
    }
    
    
    /**
     * loads subjects resource file
     */
    public static def loadSubjects(args, csvFile=null) {
        def reader = null
        
        if( csvFile ) {
            reader = new CSVReader( new File(csvFile).newReader(), ',' as Character, '"' as Character) 
        }
        else {
            def istr = Main.class.getResourceAsStream("/haplo_hla.csv")
            reader = new CSVReader(new InputStreamReader(istr), ',' as Character, '"' as Character)
        }
        
        List myEntries = reader.readAll()
        reader.close()
        
        def splitHaplo = {str-> (str.split(',') as List).collect{it.trim()} }
        
        myEntries.each{ toks->
            subjects[toks[0]] = new Subject(id:toks[0])
                
            hlaGenes.eachWithIndex{gene, i->
                subjects[toks[0]].alleles[gene] = toks[i+1] ? splitHaplo(toks[i+1]) : []
            }
        }
        
        
        //load subject lanes
        def istr = Main.class.getResourceAsStream('/'+GEUV_FILE)
        reader = new BufferedReader(new InputStreamReader(istr))
        
//        reader.readLine()//skip first line
//        reader.splitEachLine("\\s"){ toks->
//            assert toks[1]==toks[2], "Error: different lanes for ${toks[0]}"
//            
//            if( subjects.containsKey(toks[0]) ){
//                subjects[toks[0]].lanes << toks[1]
//            }
//            
//            // add to the map of all subjects
//            def newSbj = new Subject(id:toks[0])
//            newSbj.lanes << toks[1]
//            allSubjects[toks[0]] = newSbj
//        }
        
        def groups = getOptGroups(args)
        loadGeuvadis(reader, (groups ?: ['CEU']) )
        reader.close()
    }
    
    private static loadGeuvadis(reader, List groups){
        int GROUP_IDX = 6
        int FASTQ_IDX = 28
        int ID_IDX = 0
        reader.readLine()//skip first line
        
        reader.splitEachLine("\t"){ toks->
            if( toks[GROUP_IDX] in groups ){
                def id = toks[ID_IDX]
                String lane = toks[FASTQ_IDX]
                lane = lane.substring( lane.lastIndexOf('/')+1, lane.lastIndexOf('_'))

                if( subjects.containsKey(id) ){
                    if( !(lane in subjects[id].lanes) ){
                        subjects[id].lanes << lane
                    }
                }

                // add to the map of all subjects
                if( !allSubjects.containsKey(id) ){
                    def newSbj = new Subject(id:id)
                    newSbj.lanes << lane
                    allSubjects[id] = newSbj
                }
            }
        }
    }
    
    
    /**
     *
     * Generates allele sequences for each subject, the result files are named:
     * <outdir>/<id>.fa
     * Also, for each fasta file, generates the bowtie index and the fasta index.
     */
    static def generateAlleleSeq(outdir) {
        subjects.each{ id, subj-> 
            if( !id.startsWith('TEST')) {
                String sbjOutDir = outdir+id+'/'
                Utils.createDir(sbjOutDir)
                generateAlleleSeq(subj, sbjOutDir)
            } 
        }
    }
    
    static def generateAlleleSeq(subj, outdir){
        def seq1 = new StringBuilder()
        def seq2 = new StringBuilder()
        def annot1 = new StringWriter()
        def annot2 = new StringWriter()
        def prevGene = ''
        def seqIds = [subj.id+"_1", subj.id+"_2"]
        boolean reverseIfNeg = false //reverse complemented if gene is in negative strand

        hlaGenes.eachWithIndex{ gene, i->
            if( subj.alleles[gene] ) {
                //add intergenic region
                [seq1, seq2].each{ it <<  HlaSequences.getIntergenic(prevGene, gene) }
                prevGene = gene
                
                // generate annotation
                HlaSequences.writeAnnotation( gene + subj.alleles[gene][0],
                    seq1.length(), seqIds[0], annot1 )
                HlaSequences.writeAnnotation( gene + subj.alleles[gene][(subj.alleles[gene].size()==2 ? 1 : 0)],
                    seq2.length(), seqIds[1], annot2 )
                
                // get allele gene sequence
                seq1 << HlaSequences.getAlleleSeq(gene + subj.alleles[gene][0], reverseIfNeg)
                seq2 << HlaSequences.getAlleleSeq(gene + subj.alleles[gene][(subj.alleles[gene].size()==2 ? 1 : 0)], reverseIfNeg)       
                
            }
        }

        // add final intergenic region
        [seq1, seq2].each{ it << HlaSequences.getIntergenic(hlaGenes[hlaGenes.size()-1],'') }

        def faFile = outdir+subj.id+'.fa'
        def gffFile = outdir+subj.id+'.gff3'
        
        // write gff file
        def writer = new File(gffFile).newWriter()
        [annot1, annot2].each{ writer.println(it.toString()) }
        writer.close()

        def seqs = [seq1, seq2].collect{ new DNASequence(it.toString()) }
        seqs.eachWithIndex{ seq, i-> seq.setAccession(new AccessionID(seqIds[i])) }

        FastaWriterHelper.writeNucleotideSequence( new File(faFile), seqs)

        
        if( !AnnotationDB.genBowtieIndex(faFile, '2') ){ // uses bowtie2
            System.err.println("Error generating bowtie2 index for ${faFile}")
            System.exit(1)
        }
        
        if( !AnnotationDB.genFastaIndex(faFile) ){
            System.err.println("Error generating fasta index for ${faFile}")
            System.exit(1)
        }
            
    }
    
    
    /**
     *
     */
    static def hlaMap(args) {
        String outdir = args.find{ it.startsWith(OPT_OUTPUT_DIR) }
        String inputdir = args.find{ it.startsWith(OPT_INPUT_DIR) }
        def dirs = checkAndCleanDirOpts([outdir,inputdir])
        outdir = dirs[0]
        inputdir = dirs[1] //fastx_toolkit fastq files base directory
        
        def subjFile = readOption(args, OPT_SUBJECTS)
        
        //create annotation workdir
        String annotDir = outdir+'uniprotkb/'
        Utils.createDir(annotDir)
        HlaSequences.workDir = annotDir
        
        if( subjFile ) {
            // load subjects
            subjects = [:] as TreeMap
            loadSubjects(args, subjFile)
        }
        
        // load data and sequences
        HlaSequences.load()
        
        //generate hla fragments for each subject (fasta, gff3, bowtie index and fasta index)
        generateAlleleSeq(outdir)
        
        def writeHeader = { writer-> //print qsub header
            writer.println("#!/bin/bash")
            writer.println("#")
            writer.println("#PBS -N hlaMap_")
            writer.println("#PBS -M <email>")
            writer.println("#PBS -m abe")
            writer.println("#PBS -l nodes=1:ppn=7,mem=8gb")
            writer.println('#')    
        }  //      
        
        // generate tophat + cufflinks commands
        def writer = new PrintWriter( outdir+'hlamap_commands.txt')
        
        writeHeader(writer)
        
        subjects.each{ id, subj->
            if( !id.startsWith('TEST')) {
                writer.println("### Jobs for subject ${id}:")
                writer.println( JobGenerator.tophatHLA(subj, inputdir, outdir, true) )
                writer.println( JobGenerator.cufflinksHLA(subj, outdir, true) )
                writer.println('#')
            }
        }
        
        writer.close()
    }
    
    /**
     * PHLAT HLA typing
     */
    static def phlatTyping(args) {
        String outdir = args.find{ it.startsWith(OPT_OUTPUT_DIR) }
        def dirs = checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        String inputdir = readOption(args, OPT_INPUT_DIR)//fastx_toolkit fastq files base directory
        if( !inputdir.endsWith('/') ){ inputdir += '/'}
        
        def subjFile = readOption(args, OPT_SUBJECTS)
        
        def subjs = [] as TreeSet
        
        //read subjects
        def reader = new File(subjFile).newReader()
        reader.readLine()//skip header
        reader.splitEachLine("\t"){ toks->
            subjs << toks[0]
        }
        reader.close()
        
        
        def writeHeader = { writer->//print qsub header
            writer.println("#!/bin/bash")
            writer.println("#")
            writer.println("#PBS -N phlat_")
            writer.println("#PBS -M <email>")
            writer.println("#PBS -m abe")
            writer.println("#PBS -l nodes=1:ppn=7,mem=8gb")
            writer.println('#')    
        }//
        
        // generate phlat commands
        def writer = new PrintWriter( outdir+'phlat_commands.sh.txt')
        writeHeader(writer)
        
        subjs.each{ id->
            def job = 
                JobGenerator.phlatHLATyping(id, inputdir, config.org.webngs.tools.dir, config.org.webngs.scripts.dir)
            writer.println("#Job: ${job.jobId}")
            writer.println("mkdir ${StorageManager.replacePrefixs(job.url)}")
            writer.println(StorageManager.replacePrefixs(job.command))
        }
        
        writer.close()
    }
    
    /**
     * load ngsengine properties file
     */
    protected static def loadConfig(args){
        String props = readOption(args, OPT_PROPS)
        
        if( !props ){
            System.err.println("Config properties file not found")
            System.exit(1)
        }
        
        //load ngsengine properties
        config = new ConfigSlurper().parse( org.awsjoblauncher.Main.loadProperties(props) )
        
        
        StorageManager.init(config.toProperties())
        PostBase.setProps(config.toProperties())
    }
    
    /**
     * returns a closure that stores the finished jobs in a set
     */
    private static def storeFinishedJobs(){
       return {res->
            if( res.type==JobResponse.TYPE_TOBEEXEC ){
                println "[${res.jobId}] ${res.type} to be executed in ${res.machine} at ${new Date()}"
            }
            else if( res.type==JobResponse.TYPE_WASEXEC ){
                println "[${res.jobId}] ${res.type} was executed in ${res.machine} with status ${res.status} at ${new Date()}"
                finishedJobs << res
            }
       } 
    }
    
    /**
     *
     */
    protected static def execJobs(jobList){
        if( config.org.webngs.process=='local' ) {
            jobList.each{
                def res = LocalJobSlave.executeJob(it)
                println "[${it.jobId}] ${it.jobType} was executed with status ${((Integer)res[0]) ? 'Error' : 'OK'} at ${new Date()}"
                
                if( (Integer)res[0]!=0 ) { //print error message
                    System.err.println(res[1])
                }
                
                finishedJobs << res
            }
        }
        else {
            sqsMaster = QueueUtils.launchMasterQueue(config.toProperties(), storeFinishedJobs() )
            sqsReceiver = QueueUtils.launchReceiverQueue(config.toProperties())

            jobList.each{ sqsMaster.schedule(it) }

            while( true ){
                try{ Thread.sleep(5000) } catch(e){}
                if( finishedJobs.size()==jobList.size() )
                    break
            }
        }
        
        println "Exec. finished. ${finishedJobs.size()} jobs processed."
    }
    
    /**
     *
     */
    protected static def getOptSubjects(args){
        def subs = readOption(args, OPT_SUBJECTS)
        
        if( !subs )
            return null
            
        subs = subs.split(',') as List
        
        if( !subs ){
            System.err.println("Error: no subjects ids present")
            System.exit(1)
        }
        
        def selected = [:]
        
        subs.each{ 
            if( allSubjects.containsKey(it) )
                selected[it] = allSubjects[it]
        }
        
        return selected
    }
    
    /**
     *
     */
    public static def getOptGroups(args){
        def grps = readOption(args, OPT_GROUP)
        
        if( !grps )
            return null
            
        grps = grps.split(',') as List
        
        return (grps) ?: null
    }
    
    /**
     *
     */
    protected static int getStartPhase(args){
        def start = readOption(args, OPT_START_PHASE)
        
        if( !start )
            return 0
            
        try{ start = start as Integer }
        catch(e){ start = 0 }
            
        return start
    }
    
    /**
     *
     */
    static def launchQC(args){
        println "QC of fastq data in ${config.hlaqtl.ceudata.bucket}"
        def subjects = getOptSubjects(args) ?: allSubjects
        def jobList = subjects.collect{ JobGenerator.fastQC(it.value, config.hlaqtl.ceudata.bucket) }
        execJobs(jobList)
    }
    
    /**
     *
     */
    static def launchFiltering(args) {
        String inputUrl = config.hlaqtl.ceudata.bucket
        if( !inputUrl ) {
            inputUrl = StorageManager.UPLOADS_PRE
        }
        
        println "Filtering of fastq data in ${inputUrl}"
        def optMap = [minLen:'50', maxLen:'no', trimFirst:'14', qualFilt:'20', percQual:'66']
        println "Filtering parameters:\n${optMap}"
        
        def subjects = getOptSubjects(args) ?: allSubjects
        
        
        def jobList = subjects.collect{ 
                JobGenerator.fastxToolkit(it.value, inputUrl, optMap) 
            }
        execJobs(jobList)
    }
    
    /**
     *
     */
    private static def runTophat(subjects, bwtIndex, reference){
        println "... Running tophat ... ${new Date()}"

        def jobList = subjects.collect{ 
            def input = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_FASTXTOOL}_${it.value.id}/"
            JobGenerator.tophat(it.value, input, bwtIndex)
        }

        execJobs(jobList)
        finishedJobs.clear()
        println "... Finished tophat ... ${new Date()}"
    }
    
    /**
     *
     */
    private static def runCufflinks(subjects, bwtIndex, reference, boolean refOnly = false){
        if(reference){
            if( refOnly )
                println "... Running cufflinks using reference only ..."
            else
                println "... Running cufflinks RABT ..."
        }
        else{
            println "... Running cufflinks ..."
        }
        
        def jobList = subjects.collect{
            def input = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_TOPHAT}_${it.value.id}/"
            JobGenerator.cufflinks(it.value, input, reference, refOnly)
        }

        execJobs(jobList)
        finishedJobs.clear()
        println "... Finished cufflinks ... ${new Date()}"
    }
    
    /**
     * option --start:
     * (1 to 4) : 1 - tophat, 2 - cufflinks, 3 - cuffcompare, 4 - cuffcompare CEU60
     * (5 to 7) : 5 - cufflinks RABT, 6 - cuffcompare CEU60 RABT, 7 - cufflinks using reference only
     * 8 - pipeline: filtering+tophat+cufflinks (RABT, refOnly)
     */
    static def launchCBCB(args){
        def bwtIndex = config.org.webngs.bowtieIdx.dir+'hg19'
        def reference = StorageManager.ANNOTATIONS_PRE+'9606/refGene.txt.gtf.gz'
        def ensemblRef = StorageManager.ANNOTATIONS_PRE+'9606/ensGene.txt.gtf.gz'
        def subjects = getOptSubjects(args) ?: allSubjects
        int start = getStartPhase(args)
        def jobList = null
        
        if( start<=1 ){
            //Round 1 : tophat
            runTophat(subjects, bwtIndex, reference)
        }
        
        if( start<=2 ){
            //Round 2 : cufflinks
            runCufflinks(subjects, bwtIndex, null)
        }
        
        if( start<=3 ){
            //Round 3 : cuffcompare
            println "... Running cuffcompare ... ${new Date()}"

            jobList = subjects.collect{
                def input = "${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFLINKS}_${it.value.id}/"
                JobGenerator.cuffcompare(it.value, [input], reference)
            }

            execJobs(jobList)
            finishedJobs.clear()
            println "... Finished cuffcompare ... ${new Date()}"
        }
        
        if( start<=4 ){
            println "... Running cuffcompare CEU60 ... ${new Date()}"
            
            def inputList = allSubjects.collect{"${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFLINKS}_${it.value.id}/"}
            jobList = [JobGenerator.cuffcompare(new Subject(id:'CEU60'), inputList, reference)]
            execJobs(jobList)
            finishedJobs.clear()
            
            println "... Finished cuffcompare CEU60 ... ${new Date()}"
        }
        
        //RABT (reference annotation based transc. assembly)
        if( start==5 ){
            runCufflinks(subjects, bwtIndex, reference)
        }
        
        if( start==6 ){
            println "... Running cuffcompare RABT CEU60 ... ${new Date()}"
            
            def inputList = allSubjects.collect{"${StorageManager.OUTPUT_PRE}${CommConstants.COMMNAME_CUFFLINKS}_ref_${it.value.id}/"}
            jobList = [JobGenerator.cuffcompare(new Subject(id:'CEU60_ref'), inputList, null)]
            execJobs(jobList)
            finishedJobs.clear()
            
            println "... Finished cuffcompare CEU60 ... ${new Date()}"
        }
        
        // cufflinks use supplied reference
        if( start==7 ){
            runCufflinks(subjects, bwtIndex, reference, true)
        }
        
        //pipeline: tophat+cufflinks (refSeq/ensembl refOnly)
        if( start==8 ){
            runTophat(subjects, bwtIndex, reference)
            ///runCufflinks(subjects, bwtIndex, reference)
            runCufflinks(subjects, bwtIndex, reference, true)
            runCufflinks(subjects, bwtIndex, ensemblRef, true)
        }
        
        // cufflinks use ensembl reference only
        if( start==9 ){
            runCufflinks(subjects, bwtIndex, ensemblRef, true)
        }
    }
    
    
    /**
     * calling example:
     * eqtl --props=/home/victor/.grails/ngsengine.properties  
     * --output=/home/victor/Escritorio/transcriptome_60/output 
     * --locus=chr17:37124652-39366847 
     * --vcf=/home/victor/Escritorio/transcriptome_60/data/17test.vcf.gz 
     * --subjs=NA06985 --group=CEU
     */
    static def launchEqtl(args){
        //read locus, vcf and group args
        def locus = readOption(args, OPT_LOCUS)
        def vcfFile = readOption(args, OPT_VCF)
        def groups = getOptGroups(args)
        def cuffDirPrefix = readOption(args, OPT_CUFFPRE)
        def transEqtl = readOption(args, OPT_TRANSEQTLSNP)
        def isoformsFile = readOption(args, OPT_ISOFORMS)
        
        if( (!transEqtl && !locus) || !groups ){
            System.err.println( "Usage: ${USA_EQTL}")
            return null
        }
        
        if( transEqtl && !vcfFile ){
            System.err.println( "VCF dir needed!")
            return null
        }
        
        //read outdir arg
        String outdir = args.find{ it.startsWith(OPT_OUTPUT_DIR) }
        def dirs = [outdir]
        //read expresion data dir arg
        String exprdir = args.find{ it.startsWith(OPT_EXPRDATA) }
        
        if( exprdir )
            dirs << exprdir
        
        dirs = checkAndCleanDirOpts(dirs)
        outdir = dirs[0]
        
        if( exprdir )
            exprdir = dirs[1]
        
        //print readed options
        println "${LINE_PRE}Options present:"
        println "locus: ${transEqtl ? 'trans-eqtl' : locus}"
        println "groups: ${groups}"
        println "output: ${outdir}"
        println "cuffdiff pref: ${cuffDirPrefix}"
        println "vcf: ${vcfFile ?: 'will be downloaded from 1000genomes'}"
            
        if( exprdir )
            println "Expression data dir: ${exprdir}"
            
        if( transEqtl ) {
            println "SNPs file: ${transEqtl}"
            println "SNPs count: ${Utils.countLines(transEqtl)}"
        }
        
        if( isoformsFile )
            println "Isoforms file: ${isoformsFile}"
            
        //init StorageManager
        def props = getProperties()
        StorageManager.init(props)
        S3Manager.init( props.getProperty(CommConstants.AWS_ACCESS_KEY), 
            props.getProperty(CommConstants.AWS_SECRET_KEY) )
        
        //read subjects arg
        def subjects = getOptSubjects(args) ?: allSubjects
        //get results bucket from config properties
        def resBucket = config.org.webngs.aws.bucket.output
        
        // read SNPs and genotypes
        println "${LINE_PRE}Loading SNPs and genotypes"
        
        boolean _1000gOnly = true
        String snpsFile = transEqtl ?: null
        locus = transEqtl ? null : locus
        
        def snpMng = SNPManager.loadSNPData(subjects.keySet(), outdir, vcfFile,
            groups, locus, _1000gOnly, snpsFile, true)

        //read expression data of each subject
        def selIsoforms = null
        
        if( isoformsFile ){
            // keep selected isoforms only
            selIsoforms = IsoformParser.parseIsoFilterFile(isoformsFile)
        }
        
        def listIsoData = null
        File mergTrackFile = Utils.getFileIfCompress(outdir+MERG_TRACK_FILE)
        
        if( mergTrackFile==null ){
            println "${LINE_PRE}Loading expression data from ${resBucket}. ${new Date()}"
            cuffDirPrefix = (cuffDirPrefix) ?: CommConstants.COMMNAME_CUFFLINKS+'_refOnly_'
        
            listIsoData = IsoformParser.loadExpressionData( 
                subjects.values(), resBucket, cuffDirPrefix, outdir, locus, exprdir, selIsoforms)
            
            // log expression and genotypes files
            println "${LINE_PRE}Generating ${MERG_GFF_FILE} and ${MERG_TRACK_FILE} files"
            IsoformParser.writeIsoforms(listIsoData, outdir+MERG_GFF_FILE)
            IsoformParser.writeSubjectIsoforms(listIsoData, outdir+MERG_TRACK_FILE )
        }
        else{
            println "${MERG_TRACK_FILE} already present: parsing expression data"
            listIsoData = IsoformParser.parseSubjectIsoforms(mergTrackFile, subjects, selIsoforms)
        }
        
        //gather isoforms
        def isoSet = [] as TreeSet
        listIsoData.each{ 
            it.isoforms.each{key, val-> isoSet << (val.mergedId ?: val.id) } 
        }
        
        def corrCalc = new CorrelationCalc( subjects: subjects.values().collect{it.id}, 
            snps:snpMng.imputedSnps.keySet(), isoforms:isoSet, 
            threads:config.org.webngs.tool.threadCount as Integer )
        
        //populate CorrelationCalc genotypes and expression data
        subjects.values().each{ subj->
            def mapGeno = [:] as TreeMap
            def mapExpr = [:] as TreeMap
            
            //get genotypes
            snpMng.genotypes[subj.id]?.snps.each{ snpId, alleles->
                Double value = snpMng.imputedSnps[snpId].encode(alleles)
                if( value!=null )
                    mapGeno[snpId] = value
            }
            
            //get expression
            def subjIsoforms = listIsoData.find{it.subject.id==subj.id}?.isoforms
            isoSet.each{ isoId->
                if( subjIsoforms ){
                    def isoData = subjIsoforms[isoId]
                    mapExpr[isoId] = (isoData) ? isoData.fpkm : 0.0
                }
                else{
                    mapExpr[isoId] = 0.0
                }
            }
            
            corrCalc.genotypes[subj.id] = mapGeno
            corrCalc.expression[subj.id] = mapExpr
        }
        
        debugCorrData(corrCalc, outdir)
        
        println "${LINE_PRE}Perform correlation in ${corrCalc.threads} threads. ${new Date()}"
        println "${corrCalc.subjects.size()} subjects"
        println "${corrCalc.snps.size()} SNPs"
        println "${corrCalc.isoforms.size()} isoforms"
        corrCalc.calcCorrelations()
        
        //generate text results
        println "${LINE_PRE}Generating correlation results in ${CORREL_FILE}. ${new Date()}"
        corrCalc.printResults( new PrintWriter(outdir+CORREL_FILE) )
        
        // filter correlation results file
        double pvalThr = 0.05
        filterCorrResults(outdir+CORREL_FILE, MIN_OBS, pvalThr)
        
        //generate charts
        println "${LINE_PRE}Generating charts. ${new Date()}"
        def chartGen = new ChartGenerator(correlation:corrCalc, outDir:outdir, snps:snpMng.imputedSnps)
        pvalThr = 1e-3
        double corrThr = 0.4
        chartGen.generate(pvalThr, corrThr, 0)
        
        //compress results and upload to S3
        def objName
        
        if( !transEqtl ) {
            objName = compEqtlResObject(locus, cuffDirPrefix, groups)
        }
        else{
            objName = compEqtlTransObject(cuffDirPrefix, groups)
        }
        
        //println "${LINE_PRE}Uploading final results to ${resBucket}/${objName}"
        //uploadToS3(outdir, resBucket, objName, ['.png'])
        
    }
    
    /**
     * compose a result folder name for an eqtl calculation
     */
    static String compEqtlResObject(locus, cuffDirPrefix, groups) {
        def suff = locus.replace(':','_')
        suff = suff.replace('-','_')
        def cuffpart = cuffDirPrefix.substring( cuffDirPrefix.indexOf('_')+1 )
        def grppart = groups.sum{ '_'+it }
        def objName = "eqtl_${suff}_${cuffpart}${grppart}"
        
        return objName
    }
    
    static String compEqtlTransObject(cuffDirPrefix, groups) {
        def cuffpart = cuffDirPrefix.substring( cuffDirPrefix.indexOf('_')+1 )
        def grppart = groups.sum{ '_'+it }
        def date = String.format('%1$tY%1$tm%1$td',new Date())
        def objName = "eqtltrans_${date}_${cuffpart}${grppart}"
        
        return objName
    }
    
    /**
     * uploads the content of dir to S3 (not recursively)
     * 
     * @param comprExcl : list with file suffixes than wont be compressed
     */
    static protected uploadToS3(String dir, String bucket, String objBase, comprExcl = []){
        // tar folder and upload to S3
        if( dir.endsWith('/') )
            dir = dir.substring(0, dir.length()-1)
            
        def parentDir = dir.substring(0, dir.lastIndexOf('/')+1)
        def tarFileName = dir.substring(dir.lastIndexOf('/')+1)+'.tar.gz'
        
        "tar -C ${parentDir} -czf ${parentDir+tarFileName} ${dir}".execute().waitFor()
        
        def res = StorageManager.uploadToS3( bucket, parentDir+tarFileName,
                    objBase+File.separator+tarFileName )

        if( res==null )
            throw new StorageException("Error uploading to S3.");
    }
    
    /**
     * 
     * @return a CorrelationCalc object filled with the correlation result info
     */      
    static protected parseCorrFile(file){
        def doubleVal = {val->
            try{ return (val as Double) } catch(e){ return null }
        }
        
        def corrCalc = new CorrelationCalc(snps:[] as TreeSet, isoforms: [] as TreeSet,
            corrPValues:[:], corrValues:[:] as TreeMap )
        
        def reader = Utils.createReader( Utils.getFileIfCompress(file) )
        reader.readLine()
                
        reader.splitEachLine("\\s"){ toks->
            corrCalc.snps << toks[0]
            corrCalc.isoforms << toks[1]
            def lbl = corrCalc.buildLabel(toks[0],toks[1])
            corrCalc.corrValues[lbl] = doubleVal(toks[CORR_FILE_CORR])
            corrCalc.pValues[lbl] = doubleVal(toks[CORR_FILE_PVAL])
            corrCalc.corrPValues[lbl] = doubleVal(toks[CORR_FILE_CPVAL])
            corrCalc.numObs[lbl] = (toks[CORR_FILE_NOBS] as Integer)
        }
                
        reader.close()
        return corrCalc
    }
    
    /**
     * discards lines in corr. results file containing NaNs, nulls or num_obs < MIN
     */      
    static protected filterCorrResults(String file, int min, double pvalueThr){
        def writer = new PrintWriter(file+'.filter')
        def reader = (new File(file)).newReader()
        
        writer.println( reader.readLine() ) //write header
        
        reader.splitEachLine("\\s"){ toks->
            //correlation value
            if( !(toks[CORR_FILE_CORR] in ['NaN','null']) ){
                try{
                    def pval = toks[CORR_FILE_CPVAL] as Double
                    def nobs = toks[CORR_FILE_NOBS] as Integer
                    
                    if( pval<pvalueThr && nobs>min ){
                        writer.println(
                            "${toks[0]}\t${toks[1]}\t${toks[2]}\t${toks[3]}\t${toks[4]}\t${toks[5]}\t${toks[6]}\t${toks[7]}"
                        )
                    }
                } catch(e){}
            }
        }
        
        reader.close()
        writer.close()
    }
    
    /**
     *
     */
    static protected debugCorrData(CorrelationCalc corrData, outdir){
//        def subjects //subjects ids
//        def snps //SNPs ids
//        def isoforms //Isoforms ids
//        Map genotypes = [:] as TreeMap //map with key=subject value=map of pairs {SNP_id, value}
//        Map expression = [:] as TreeMap //map with key=subject value=map of pairs {Isoform_id, value}

        def writer = new PrintWriter(outdir+CORR_DEBUG_FILE)
        
        corrData.subjects.each{ subjId->
            corrData.genotypes[subjId]?.each{ snp, val->
                writer.println("${subjId}\t${snp}\t${val}")
            }
            corrData.expression[subjId]?.each{ iso, val->
                writer.println("${subjId}\t${iso}\t${val}")
            }
        }
        
        writer.close()
        
        "gzip -f ${outdir+CORR_DEBUG_FILE}".execute().waitFor()
    }
    
    /**
     *
     */
    static protected downExpresData(args) {
        def cuffDirPrefix = readOption(args, OPT_CUFFPRE)
        def subjects = allSubjects
        def groups = getOptGroups(args)
        
        //read outdir arg
        String outdir = args.find{ it.startsWith(OPT_OUTPUT_DIR) }
        def dirs = checkAndCleanDirOpts([outdir])
        outdir = dirs[0]
        
        if( !groups ){
            System.err.println( "Usage: ${USA_DOWN_EXPR}")
            return null
        }
        
        //get results bucket from config properties
        def resBucket = config.org.webngs.aws.bucket.output
        
        println "Download expression data from:${resBucket} prefix:${cuffDirPrefix}"
        println "groups: ${groups}"
        println "output: ${outdir}"
       
        cuffDirPrefix = (cuffDirPrefix) ?: CommConstants.COMMNAME_CUFFLINKS+'_refOnly_'
        
        IsoformParser.downloadExpressionData( subjects.values(), resBucket, cuffDirPrefix, outdir)
    }
    
    
}
