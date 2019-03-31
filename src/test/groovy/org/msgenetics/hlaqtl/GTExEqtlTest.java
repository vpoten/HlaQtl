/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import java.util.ArrayList;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import tech.tablesaw.api.Table;
import tech.tablesaw.api.DoubleColumn;

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
        assertEquals(tissues.size(), 48);
    }
    
    @Test
    public void testLoadTissue() {
        GTExEqtl instance = new GTExEqtl();
        instance.setPath(dataPath);
        Table table = instance.loadTable("Adipose_Subcutaneous");
        assertNotNull(table);
        
        assertEquals(table.rowCount(), 23935);
        
        List<String> names = table.columnNames();
        assertEquals(names.get(2), "gene_chr");
        assertEquals(names.get(14), "pos");
        assertEquals(names.get(18), "rs_id_dbSNP147_GRCh37p13");
        assertEquals(names.get(27), "pval_beta");
        assertEquals(names.get(28), "qval");
        
        DoubleColumn column = (DoubleColumn) table.column(27);
        assertTrue(column.max() < 1);
        assertTrue(column.min() > 0);
        
        Table filtered = table.where(column.isLessThan(0.05));
        int num = filtered.rowCount();
        assertTrue(num < 23935);
        assertTrue(num > 0);
    }
    
    @Test
    public void testLoadAllTissues() {
        GTExEqtl instance = new GTExEqtl();
        instance.setPath(dataPath);
        
        List<String> tissues = (List<String>) instance.getTissues();
        List<Table> tables = new ArrayList<Table>();
        
        for(String tissue : tissues) {
            Table table = instance.loadTable(tissue);
            assertNotNull(table);
            assertTrue(table.rowCount() > 1000);
            tables.add(table);
        }
        
        assertEquals(tables.size(), 48);
    }
    
    @Test
    public void testGetBestEqtlsAllTissues() {
        Table table = GTExEqtl.getBestEqtlsAllTissues(dataPath, 0.01);
        assertTrue(table.rowCount() > 1000);
        assertTrue(table.name().startsWith("best-eqtls"));
        
        Table table2 = GTExEqtl.getBestEqtlsAllTissues(dataPath, 0.001);
        assertTrue(table.rowCount() > table2.rowCount());
    }
}
