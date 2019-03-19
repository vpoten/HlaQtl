/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class SimpleBLASTTest {
    
    private static String BLAST_PATH = "/home/victor/software/ncbi-blast-2.2.28+/bin";
    
    public SimpleBLASTTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void runBlastn() {
        SimpleBLAST instance = new SimpleBLAST(BLAST_PATH);
        
        String query = "/home/victor/software/ncbi-blast-2.2.28+/fasta/urod/mouse.fasta";
        String subject = "/home/victor/software/ncbi-blast-2.2.28+/fasta/urod/human.fasta";
        
        BLASTResult res = instance.runBlastn(query, subject);
        assertNotNull(res);
        AlignResult ali = res.getAlignments().get(0);
        assertNotNull( ali.getExpect() );
        assertTrue( res.getQuery().startsWith("ENSMUSG") );
        assertTrue( res.getSubject().startsWith("ENSG") );
        Integer exp = 738;
        assertEquals( ali.getScore(), exp);
        exp = 1104;
        assertEquals( res.getSubjectLength(), exp);
        assertEquals( ali.getIdentities(), res.averageIdent());
        assertEquals( ali.getScore(), res.totalScore());
        
        subject = "/home/victor/software/ncbi-blast-2.2.28+/fasta/urod/zebrafish.fasta";
        res = instance.runBlastn(query, subject);
        assertNotNull(res);
        assertTrue( res.getAlignments().isEmpty() );
    }
    
}