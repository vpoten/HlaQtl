/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import com.mashape.unirest.http.exceptions.UnirestException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import org.ngsutils.variation.SNPData;

/**
 *
 * @author victor
 */
public class EnsemblRestClientTest {
    
    public EnsemblRestClientTest() {
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
    public void getSnps() throws UnirestException, InterruptedException {
        ArrayList<String> ids = new ArrayList<String>();
        ids.add("rs12722489");
        ids.add("rs6897932");
        ids.add("rs6498169");
        ids.add("rs6604026");
        ids.add("rs10984447");
        
        EnsemblRestClient instance = new EnsemblRestClient("http://grch37.rest.ensembl.org", 15, 200);
        List<SNPData> result = instance.getSnps(ids, "human");
        assertNotNull(result);
        assertEquals(result.size(), ids.size());
        
        HashSet<String> idsSet = new HashSet<String>(ids);
        
        for(int i=0; i<result.size(); i++) {
            assertTrue(idsSet.contains(result.get(i).getId()));
        }
    }
}
