/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import tech.tablesaw.api.Table;

/**
 *
 * @author victor
 */
public class SimplifiedEqtlTest {
    
    public SimplifiedEqtlTest() {
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
    public void testLoadTable() {
        SimplifiedEqtl instance = new SimplifiedEqtl();
        instance.setPath("/home/victor/Escritorio/transcriptome_60/eqtl_search/eqtl_simplified_sample.csv");
        Table table = instance.loadTable("eqtls");
        assertTrue(table.rowCount() > 1900);
        
        table = instance.getBestEqtls(1e-4);
        assertTrue(table.rowCount() < 1400);
    }
}
