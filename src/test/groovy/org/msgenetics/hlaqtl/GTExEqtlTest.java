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

import tech.tablesaw.api.Table;

/**
 *
 * @author victor
 */
public class GTExEqtlTest {
    
    static String dataPath = "/home/victor/Descargas/GTEx_Analysis_v7_eQTL";
    
    public GTExEqtlTest() {
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
    public void testGetTissues() {
        GTExEqtl instance = new GTExEqtl();
        instance.setPath(dataPath);
        List<String> tissues = (List<String>) instance.getTissues();
        assertEquals(tissues.size(), 46);
    }
    
    @Test
    public void testLoadTissue() {
        GTExEqtl instance = new GTExEqtl();
        instance.setPath(dataPath);
        Table table = instance.loadTable("Adipose_Subcutaneous");
        assertNotNull(table);
    }
}
