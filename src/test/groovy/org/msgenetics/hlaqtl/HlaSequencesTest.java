/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.msgenetics.hlaqtl.hla.HLAAllele;

/**
 *
 * @author victor
 */
public class HlaSequencesTest {
    
    static boolean TEST_LOAD_GENE = true; //5 minute tests (first time is called)
    
    public HlaSequencesTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        HlaSequences.setWorkDir("/home/victor/Escritorio/tests_ngsengine/uniprotkb/");
        HlaSequences.loadPositions();
        HlaSequences.loadIntergenic();
        
        
        if( TEST_LOAD_GENE ){
            HlaSequences.loadGene((String) Main.getDRB1());
            HlaSequences.loadGene((String) Main.getDRB5());
            HlaSequences.loadGene((String) Main.getDQB1());
        }
        
        if( Main.getSubjects().isEmpty() ) {
            Main.loadSubjects( new ArrayList<String>() );
        }
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void downloadDna() {
        System.out.println("TestDownloadDna");
        
        Map seqs = (Map) HlaSequences.downloadDna("chr1", 1000000, 1001000, "hg19");
        assertEquals(seqs.size(), 1);
        
        String key = "hg19_dna range=chr1:1000001-1001000 5'pad=0 3'pad=0 strand=+ repeatMasking=none";
        assertTrue( seqs.containsKey(key) );
        
        DNASequence seq = (DNASequence) seqs.get(key);
        assertTrue( seq.getSequenceAsString().endsWith("GAGAGATTGTAATAAATAAAGAC") );
        assertEquals( seq.getSequenceAsString().length(), 1000);
    }
    
    @Test
    public void haplotypeSet() {
        System.out.println("TestHaplotypeSet");
        
        Set haplo = HlaSequences.haplotypeSet(Main.getSubjects().values());
        assertTrue( !haplo.isEmpty() );
        
        Set geno = HlaSequences.genotypeSet(Main.getSubjects().values());
        assertTrue( !geno.isEmpty() );
        
        Map map = HlaSequences.countGenotypes(Main.getSubjects().values());
        assertTrue( !map.isEmpty() );
        
    }
    
    
    @Test
    public void load() {
        System.out.println("TestLoad");
        
        Map posMap = (Map) HlaSequences.getPositionMap();
        Map geneMap = (Map) posMap.get(Main.getDRB1());
        
        assertEquals( geneMap.get("position"), "chr6:32546547-32557613" );
        assertEquals( geneMap.get("strand"), "-" );
        
        
    }
    
    @Test
    public void getIntergenic() {
        System.out.println("TestGetIntergenic");
        
        DNASequence seq = (DNASequence) HlaSequences.getIntergenic( Main.getDRB5(), Main.getDRB1());
        assertNotNull(seq);
        
        seq = (DNASequence) HlaSequences.getIntergenic( Main.getDQB1(), "");
        assertNotNull(seq);
        
    }
    
    @Test
    public void getAllele(){
        if( TEST_LOAD_GENE ) {
            System.out.println("TestGetAllele");

            String name = HlaSequences.getAlleleName("DRB5*01:01:01");
            assertEquals(name, "DRB5*01:01:01");

            name = HlaSequences.getAlleleName("DRB1*09:30");
            assertEquals(name, "DRB1*09:01:02");

            DNASequence seq = (DNASequence) HlaSequences.getAlleleSeq("DRB5*01:01:01", true);
            assertNotNull(seq);

            DNASequence seq2 = (DNASequence) HlaSequences.getAlleleSeq("DRB5*01:01:01", false);

            String str1 = seq.toString();
            String str2 = seq2.toString();

            assertTrue( !str1.equals(str2) );

            seq = (DNASequence) HlaSequences.getAlleleSeq("DRB1*03:01:11", true);
            assertNotNull(seq);

            //check HLAAllele objects
            HLAAllele allObj = HlaSequences.getAllele("DRB5*02:06");
            assertEquals(allObj.getHlaId(),"HLA08085");
            assertEquals(allObj.getName(),"DRB5*02:06");
            assertTrue(allObj.getIsPartial());
            assertEquals(allObj.getParent().getName(),"DRB5*02:02");
            assertNotNull( allObj.getGeneSeq() );


            allObj = HlaSequences.getAllele("DRB1*03:01:03");
            assertEquals(allObj.getHlaId(),"HLA02769");
            assertEquals(allObj.getName(),"DRB1*03:01:03");
            assertTrue(allObj.getIsPartial());
            assertEquals(allObj.getParent().getName(),"DRB1*03:01:01:01");
            assertNotNull( allObj.getGeneSeq() );

            allObj = HlaSequences.getAllele("DRB1*03:01:01:01");
            assertEquals(allObj.getHlaId(),"HLA00671");
            assertEquals(allObj.getName(),"DRB1*03:01:01:01");
            assertFalse(allObj.getIsPartial());
            assertNull( allObj.getParent() );


            // check annotation
            List exons = HlaSequences.getAnnotMap().get(Main.getDRB5()).getIsoforms().get("DRB5*02:06").getExons();

            assertEquals(exons.size(), 6);

            //check sequences
            allObj = HlaSequences.getAllele("DRB1*01:54");
            assertTrue(allObj.getIsPartial());
            assertEquals( allObj.getCds().getLength(), 270);
            assertEquals( allObj.getParent().getGeneSeq().getLength(), allObj.getGeneSeq().getLength());
            assertEquals( allObj.getGeneSeq().getLength(), 10741);
            assertFalse( allObj.getParent().getGeneSeq().getSequenceAsString().equals(
                allObj.getGeneSeq().getSequenceAsString()) );

            allObj = HlaSequences.getAllele("DRB5*01:08N");
            assertTrue(allObj.getIsPartial());
            assertEquals( allObj.getCds().getLength(), 500);
            assertTrue( allObj.getParent().getGeneSeq().getLength() > allObj.getGeneSeq().getLength());
            assertTrue( allObj.getGeneSeq().getLength()>1000);
        }
        else {
            System.out.println("Skipping TestGetAllele");
        }
    }
    
    @Test
    public void generateAlleleSeq() throws IOException, CompoundNotFoundException {
        String outdir = "/home/victor/Escritorio/tests_ngsengine/";
        
        if( TEST_LOAD_GENE ) {
            System.out.println("TestGenerateAlleleSeq");
            
            Map subjects = (Map) Main.getSubjects();
            Main.generateAlleleSeq( subjects.get("TESTSBJ"), outdir);

            assertTrue( new File(outdir+"TESTSBJ.fa").exists() );
            assertTrue( new File(outdir+"TESTSBJ.gff3").exists() );
            
            Reader reader = org.ngsutils.Utils.createReader(new File(outdir+"TESTSBJ.fa"));
            
            // load secuences
            Map<String,DNASequence> sequences = (Map<String,DNASequence>) HlaSequences.readFasta(reader);
            // load annotations
            FeatureList ftList = GFF3Reader.read(outdir+"TESTSBJ.gff3");
            
            
            HLAAllele bestAllele = 
                    getBestAllele(sequences.get("TESTSBJ_1"), ftList,"TESTSBJ_1", (String) Main.getDRB1(), "DRB1*04:01:01");
            ///HlaSequences.printSequence(Main.getDRB1(), "DRB1*04:01:01");
            
            assertEquals( bestAllele.getName(), "DRB1*04:01:01");
            
            
            bestAllele = 
                    getBestAllele(sequences.get("TESTSBJ_1"), ftList,"TESTSBJ_1", (String) Main.getDRB5(), "DRB5*01:01:01");
            
            assertEquals( bestAllele.getName(), "DRB5*01:01:01");
            
            bestAllele = 
                    getBestAllele(sequences.get("TESTSBJ_1"), ftList,"TESTSBJ_1", (String) Main.getDQB1(), "DQB1*02:01:01");
            
            assertEquals( bestAllele.getName(), "DQB1*02:01:01");
            
            bestAllele = 
                    getBestAllele(sequences.get("TESTSBJ_2"), ftList,"TESTSBJ_2", (String) Main.getDRB1(), "DRB1*15:01:01:01");
            
            assertEquals( bestAllele.getName(), "DRB1*15:01:01:01");
            
            bestAllele = 
                    getBestAllele(sequences.get("TESTSBJ_2"), ftList,"TESTSBJ_2", (String) Main.getDRB5(), "DRB5*01:13");
            
            assertEquals( bestAllele.getName(), "DRB5*01:13");
        
        }
        else {
            System.out.println("Skipping TestGenerateAlleleSeq");
        }
        
    }
    
    /**
     * 
     * @param seq
     * @param ftList
     * @param seqName
     * @param gene
     * @param alleName
     * @return
     * @throws CompoundNotFoundException 
     */
    protected HLAAllele getBestAllele(DNASequence seq, FeatureList ftList, String seqName, String gene, String alleName) 
            throws CompoundNotFoundException {
        String codedName = alleName.replace(":", "%3A")+'_'+seqName;
        FeatureList subList = ftList.selectByType("CDS").selectByAttribute("Parent", codedName);
            
        StringBuilder strb = new StringBuilder();

        for(FeatureI ft : subList) { //compose CDS from features
            if( ft.seqname().equals(seqName) ){
                strb.append(  seq.getSequenceAsString(
                            ft.location().bioStart(), ft.location().bioEnd(), Strand.POSITIVE) 
                );
            }
        }

        Sequence target = new DNASequence(strb.toString(), DNACompoundSet.getDNACompoundSet());

        ///ArrayList<String> included = new ArrayList<String>();
        ///included.add(alleName);
        Map aligners = (Map) HlaSequences.alignAlleles(gene, target, 0, null);
        
        return HlaSequences.getBestAllele(aligners);
    }
    
}
