/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author victor
 */
public class ImmunobaseAssocTest {
    
    public ImmunobaseAssocTest() {
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
    
//    @Test
//    public void regionMatch() throws FileNotFoundException {
//        String tabFile = "/home/victor/Escritorio/transcriptome_60/immunochip/IC_RESULTS_15_02_2013.tab.gz";
//        String eqtlFile = "/home/victor/Escritorio/eqtls_snps/all_eqtl_cis_trans.filtered_05.txt.gz";
//        String snpsFile = "/home/victor/Escritorio/eqtls_snps/snps_info_05.log";
//        String isoFile = "/home/victor/Escritorio/eqtls_snps/filtered_isoforms_05.final.txt";
//        
//        
//        ImmunobaseAssoc instance = new ImmunobaseAssoc();
//        instance.setUsePval(false);
//        instance.setIsoGeneMap( IsoformParser.parseIsoFilterFileGenes(isoFile) );
//        instance.parseSnpsLog(snpsFile);
//        instance.parseEqtlFile(eqtlFile);
//        
//        instance.parseImmunobaseTab(tabFile);
//        Map res = instance.regionMatch("chr2:60868614-62057317"/*,"/home/victor/Escritorio/"*/);
//        assertNotNull(res);
//    }
    
//    @Test
//    public void regionGff3Match() throws FileNotFoundException {
//        boolean selectByLD = true;
//        boolean genLDgroups = false;
//        
//        String tabFile = "/home/victor/Escritorio/transcriptome_60/immunochip/IC_RESULTS_15_02_2013.tab.gz";
//        String eqtlFile = "/home/victor/Escritorio/eqtls_snps/all_eqtls_cis_nonsig.filtered.txt.gz";
//        String snpsFile = "/home/victor/Escritorio/eqtls_snps/snps_info_05.log";
//        String isoFile = "/home/victor/Escritorio/eqtls_snps/filtered_isoforms_05.final.txt";
//        
//        // MS regions data
//        String gff3Dir = "/home/victor/Escritorio/transcriptome_60/immunochip/MS/";
//        
//        String [] loci = new String [] {
//            "chr12:123316881-124008073", "chr1:2363327-2804473",
//            "chr12:57867726-58534768", "chr12:9436282-9973552",
//            "chr14:88217284-88649410", "chr16:11017058-11466511",
//            "chr2:231045388-231276234", "chr8:79161910-79760455"
//        };
//        
//        ImmunobaseAssoc instance = new ImmunobaseAssoc();
//        instance.setUsePval(false);
//        instance.setIsoGeneMap( IsoformParser.parseIsoFilterFileGenes(isoFile) );
//        instance.parseSnpsLog(snpsFile);
//        
//        //get snps set from gff3 files
//        TreeMap<String,Double> snps = new TreeMap<String,Double>();
//        for(String locus : loci) {
//            String file = gff3Dir+locus.replace(":","_").replace("-","_")+".gff3";
//            instance.parseGff3(file, snps);
//        }
//        
//        instance.parseEqtlFileBySnps(eqtlFile, 1.0, snps.keySet());
//        
//        if( genLDgroups ){
//            for(String locus : loci) {
//                String file = gff3Dir+locus.replace(":","_").replace("-","_")+".gff3";
//                instance.addGff3ToLDGroups(file, locus, gff3Dir+"ld_test/");
//            }
//            
//            instance.writeLdSnpsGroups(gff3Dir+"ms_ld_groups.txt");
//        }
//        else{
//            HashMap<Integer, PrintWriter> writerMap = new HashMap<Integer, PrintWriter>();
//            writerMap.put(ImmunobaseAssoc.getDIST_EU(), new PrintWriter(gff3Dir+"assoc.eucl.result") );
//            writerMap.put(ImmunobaseAssoc.getDIST_PW(), new PrintWriter(gff3Dir+"assoc.pwr.result") );
//            writerMap.put(ImmunobaseAssoc.getDIST_LW(), new PrintWriter(gff3Dir+"assoc.lin.result") );
//
//            for(PrintWriter writer : writerMap.values()){
//                instance.printResultHeader(writer);
//            }
//
//            for(String locus : loci) {
//                String file = gff3Dir+locus.replace(":","_").replace("-","_")+".gff3";
//                Map res = null;
//
//                if(selectByLD){
//                    res = instance.gff3MatchLD(file, "MS", locus, gff3Dir+"ld_test/");
//                }
//                else{
//                    res = instance.gff3Match(file, "MS");
//                }
//
//                for(Integer distType : writerMap.keySet()){
//                    instance.printResult(res, writerMap.get(distType), locus, 25, distType);
//                }
//            }
//
//            for(PrintWriter writer : writerMap.values()){
//                writer.close();
//            }
//        }
//    }
    
}