/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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
    public void perform() {
        String [] populations = {"CEU"};
        // String [] chrs = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
        // String [] chrs = {"X"};
        
        GTExSearcher instance = new GTExSearcher();
        instance.setSubjects(populations);
        instance.setWorkDir(workDir);
        instance.setGenomesDir(genomesDir);
        instance.setGtexDir(gtexDir);
        // instance.setChrAllowed(chrs);
        instance.loadQuerySnpsFromFile(snpsFile);
        instance.setUseCache(false);
        
        assertTrue(instance.getQueryIds().size() > 100);
        assertTrue(instance.getSubjects().size() > 100);
        
        instance.perform();
        
        // TODO
    }
}
