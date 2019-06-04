/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import java.util.ArrayList;
import java.util.Arrays;
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
public class GTExSearcherTest {
    static String workDir;
    static String gtexDir;
    static String genomesDir;
    static String snpsFile;
    
    public GTExSearcherTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        GTExSearcherTest.workDir = "/home/victor/Escritorio/gtex_eqtl_workdir";
        GTExSearcherTest.gtexDir = "/home/victor/Descargas/GTEx_Analysis_v7_eQTL";
        GTExSearcherTest.genomesDir = "/home/victor/1000genomes";
        GTExSearcherTest.snpsFile = "/home/victor/Escritorio/gtex_eqtl_workdir/MS.txt";
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
    public void createFromArgs() {
        String [] commandParts = {
            "--snps=" + GTExSearcherTest.snpsFile,
            "--workDir=" + GTExSearcherTest.workDir,
            "--gtexDir=" + GTExSearcherTest.gtexDir,
            "--genomesDir=" + GTExSearcherTest.genomesDir,
        };
        
        for (String c : commandParts)  {
            // check error with missing parameters
            String [] args = {c};
            GTExSearcher instance = GTExSearcher.createFromArgs(args);
            assertNull(instance);
        }
        
        // check a valid call with defaults
        GTExSearcher instance = GTExSearcher.createFromArgs(commandParts);
        assertNotNull(instance);
        assertNotNull(instance.getWorkDir());
        assertNotNull(instance.getGtexDir());
        assertNotNull(instance.getGenomesDir());
        assertTrue(instance.getQueryIds().size() > 100);
        assertTrue(instance.getSubjects().size() > 100);
        assertEquals(instance.getSnpRegionSize(), 10000000);
        assertTrue(Math.abs(0.05d - instance.getEqtlThr()) < 1e-12);
        assertTrue(Math.abs(0.01d - instance.getMaf()) < 1e-12);
        
        // check an invalid call
        String [] args = Arrays.copyOf(commandParts, commandParts.length + 1);
        args[commandParts.length] = "--regionSize=axx";
        instance = GTExSearcher.createFromArgs(args);
        assertNull(instance);
        
        // check a call with non-defaults
        args = Arrays.copyOf(commandParts, commandParts.length + 4);
        args[commandParts.length] = "--eqtlThr=0.01";
        args[commandParts.length + 1] = "--ldThr=0.7";
        args[commandParts.length + 2] = "--regionSize=500000";
        args[commandParts.length + 3] = "--maf=0.02";
        instance = GTExSearcher.createFromArgs(args);
        assertTrue(Math.abs(0.01d - instance.getEqtlThr()) < 1e-12);
        assertTrue(Math.abs(0.7d - instance.getLdThr()) < 1e-12);
        assertTrue(Math.abs(0.02d - instance.getMaf()) < 1e-12);
        assertEquals(instance.getSnpRegionSize(), 500000);
    }
    
//    @Test
//    public void perform() {
//        String [] populations = {"CEU"};
//        // String [] chrs = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
//        // String [] chrs = {"X"};
//        
//        GTExSearcher instance = new GTExSearcher();
//        instance.setSubjects(populations);
//        instance.setWorkDir(workDir);
//        instance.setGenomesDir(genomesDir);
//        instance.setGtexDir(gtexDir);
//        // instance.setChrAllowed(chrs);
//        instance.loadQuerySnpsFromFile(snpsFile);
//        instance.setUseCache(false);
//        
//        assertTrue(instance.getQueryIds().size() > 100);
//        assertTrue(instance.getSubjects().size() > 100);
//        
//        instance.perform();
//        
//        // TODO
//    }
}
