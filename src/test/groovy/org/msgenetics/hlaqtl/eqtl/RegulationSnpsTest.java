/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import org.msgenetics.hlaqtl.Main;
import java.util.ArrayList;
import java.util.List;
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
public class RegulationSnpsTest {
    
    public RegulationSnpsTest() {
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
    
//    @Test
//    public void generatePatterns() {
//        String workDir = "/home/victor/Escritorio/transcriptome_60/regulation_snps/";
//        String snpsInfoFile = "/home/victor/Escritorio/eqtls_snps/snps_info_05.log";
//        String cellLine = "Gm12878";
//                
//        List patterns = RegulationSnps.generatePatterns(workDir, snpsInfoFile, cellLine);
//        RegulationSnps.writePatterns(workDir+"snps_regulation.txt", patterns);
//    }
    
//    @Test
//    public void groupRanking() {
//        ArrayList<String> args = new ArrayList<String>();
//        args.add(Main.OPT_OUTPUT_DIR+"/home/victor/Escritorio/transcriptome_60/regulation_snps/");
//        args.add(Main.OPT_SNPS+"/home/victor/Escritorio/transcriptome_60/regulation_snps/snps_regulation.raw.txt.gz");
//        args.add(Main.OPT_GROUP+"/home/victor/Escritorio/eqtls_snps/group_besteqtls_0.05.txt.gz");
//        
//        RegulationSnps.performGroupRanking(args);
//    }
    
    
//    @Test
//    public void transfacEnrichment() {
//        ArrayList<String> args = new ArrayList<String>();
//        args.add(Main.OPT_OUTPUT_DIR+"/home/victor/Escritorio/transcriptome_60/regulation_snps/");
//        args.add(Main.OPT_SNPS+"/home/victor/Escritorio/transcriptome_60/regulation_snps/snps_regulation.raw.txt.gz");
//        args.add(Main.OPT_VALUE+"/home/victor/Escritorio/eqtls_snps/group_besteqtls_0.05_snpsOnly.txt");
//        
//        RegulationSnps.performTransfacEnrichment(args);
//    }
    
    
}
