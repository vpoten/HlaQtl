/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import java.util.ArrayList;
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
public class SnpGeneRelationTest {
    
    public SnpGeneRelationTest() {
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
    public void perform() {
        ArrayList<String> args = new ArrayList<String>();
        
        args.add(Main.OPT_OUTPUT_DIR+"/home/victor/Escritorio/transcriptome_60/snp_gene_out");
        args.add(Main.OPT_INPUT_DIR+"/home/victor/Escritorio/transcriptome_60");
        args.add(Main.OPT_TAXID+"9606");
        args.add(Main.OPT_SNPS+"/home/victor/Escritorio/transcriptome_60/Lista_MS.txt");
        args.add(Main.OPT_PAIRS+"/home/victor/Escritorio/transcriptome_60/NIHMS36028_snp_gene.txt");
        
        SnpGeneRelation.perform(args);
    }
    
}