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
    
    public GTExSearcherTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        // TODO
        GTExSearcherTest.workDir = "";
        GTExSearcherTest.gtexDir = "";
        GTExSearcherTest.genomesDir = "";
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
        ArrayList<String> snpIds = new ArrayList<String>();
        
        GTExSearcher instance = new GTExSearcher();
        instance.setSubjects(populations);
        instance.setWorkDir(workDir);
        instance.setGenomesDir(genomesDir);
        instance.setGtexDir(gtexDir);
        instance.setQueryIds(snpIds);
        
        instance.perform();
        
        // TODO
    }
}
