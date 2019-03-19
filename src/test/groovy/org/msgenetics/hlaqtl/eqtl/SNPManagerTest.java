/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import org.ngsutils.variation.SNPData;
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
public class SNPManagerTest {
    
    public SNPManagerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
     
    
    @Test
    public void testSNPData() {
        SNPData instance = new SNPData();
        instance.setAlleles("AT");
        
        assertTrue( Math.abs(instance.encode("AA")-0.0) < 1e-9 );
        assertTrue( Math.abs(instance.encode("AT")-1.0) < 1e-9 );
        assertTrue( Math.abs(instance.encode("TT")-2.0) < 1e-9 );
        assertNull( instance.encode("0T") );
         
        assertEquals( instance.decode(0), "AA");
        assertEquals( instance.decode(1), "AT");
        assertEquals( instance.decode(2), "TT");
                
    }
    
}
