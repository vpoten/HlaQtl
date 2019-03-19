/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.msgenetics.hlaqtl.Main;
import org.ngsutils.variation.SNPData;

/**
 *
 * @author victor
 */
public class LDCalculatorTest {
    
    public LDCalculatorTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        Main.loadSubjects( new ArrayList<String>() );
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
    
    static boolean equalDouble(double a, double b){
        return Math.abs(a-b)<1e-3;
    }

//    @Test
//    public void test() {
//        String locus = "chr12:123316881-124008073";
//        
//        ArrayList<String> snpsQuery = new ArrayList<String>();
//        snpsQuery.add("rs2343685");
//        snpsQuery.add("rs185474961");
//        
//        ArrayList<String> subjects = new ArrayList<String>();
//        subjects.addAll(Main.getAllSubjects().keySet());
//        
//        String outdir = "/home/victor/Escritorio/transcriptome_60/immunochip/MS/ld_test/";
//        
//        Map<String,Map<String,Double>> result = LDCalculator.perform(locus, snpsQuery, null, subjects, outdir);
//        assertNotNull(result);
//    }
    
    @Test
    public void test2() {
        String tpedFile = "/home/victor/Escritorio/tests_ngsengine/1000g_snps_chr1.tped";
        
        ArrayList<String> snps = new ArrayList<String>();
        String [] snp = new String [] { "rs6676197", "rs11166560", 
            "rs7533072", "rs4907906"};
        snps.addAll(Arrays.asList(snp));
       
        
        Map<String,SNPData> mapSnps = SNPData.createFromTped(tpedFile, snps);
        
        Double result = LDCalculator.calcLD(mapSnps.get(snp[1]),mapSnps.get(snp[3]));
        assertTrue( equalDouble(result,0.001) );
        
        result = LDCalculator.calcLD(mapSnps.get(snp[2]),mapSnps.get(snp[3]));
        assertTrue( equalDouble(result,0.182) );
        
        result = LDCalculator.calcLD(mapSnps.get(snp[0]),mapSnps.get(snp[1]));
        assertTrue( equalDouble(result,0.003) );
    }
    
//    @Test
//    public void test3() {
//        String tpedFile = "/home/victor/Escritorio/transcriptome_60/genotypes_motifAndGwas/1000g_snps_chr4.tped";
//        
//        ArrayList<String> snps = new ArrayList<String>();
//        String [] snp = new String [] { "rs10007052", "rs11131799", "rs11248061"};
//        snps.addAll(Arrays.asList(snp));
//        
//        Map<String,SNPData> mapSnps = SNPData.createFromTped(tpedFile, snps);
//        
//        Double result = LDCalculator.calcLD(mapSnps.get(snp[0]),mapSnps.get(snp[1]));
//        assertTrue( equalDouble(result,0.002) );
//        
//        result = LDCalculator.calcLD(mapSnps.get(snp[0]),mapSnps.get(snp[2]));
//        assertTrue( equalDouble(result,0.0001) );
//        
//    }
    
}
