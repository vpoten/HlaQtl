/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

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
public class IGSRSamplesTest {
    
    public IGSRSamplesTest() {
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
    public void create() {
        IGSRSamples instance = IGSRSamples.create();
        assertNotNull(instance);
        String [] populations = {"CEU"};
        List<String> result = instance.getSubjects(populations);
        assertNotNull(result);
        assertTrue(result.size() > 100);
    }
}
