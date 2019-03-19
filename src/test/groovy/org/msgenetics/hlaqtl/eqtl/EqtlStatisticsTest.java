/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import java.util.ArrayList;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.msgenetics.hlaqtl.Main;

/**
 *
 * @author victor
 */
public class EqtlStatisticsTest {
    
    public EqtlStatisticsTest() {
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
        
        args.add(Main.OPT_OUTPUT_DIR+"/home/victor/Escritorio/transcriptome_60/eqtl_stats_out");
        args.add(Main.OPT_INPUT_DIR+"/home/victor/Escritorio/transcriptome_60");
        args.add(Main.OPT_BLAST+"/home/victor/software/ncbi-blast-2.2.28+/bin");
        args.add(Main.OPT_SNPS+"/home/victor/Escritorio/eqtls_snps/best_eqtls_1e05.raw.out");
        args.add(Main.OPT_EQTL_DIR+"/home/victor/Escritorio/transcriptome_60");
        
        //EqtlStatistics.perform(args);
        //EqtlStatistics.performSearchCisTrans(args);
    }
    
}