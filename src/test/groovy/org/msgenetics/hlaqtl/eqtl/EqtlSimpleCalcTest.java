/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl.eqtl;

import org.msgenetics.hlaqtl.Main;
import java.util.ArrayList;
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
public class EqtlSimpleCalcTest {
    
    public EqtlSimpleCalcTest() {
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
//    public void calcEqtls() {
//        String path = "/home/victor/Escritorio/eqtl_VDR/";
//        String outdir = path+"output";
//        String vcfFile = path+"12.47770800-48770800.ALL.chr12.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz";
//        String exprFile = path+"expression_VDR.csv";
//        String locus = "chr12:47770800-48770800";
//        boolean _1000gOnly = false;
//        boolean _hapMapOnly = true;
//        
//        ArrayList<String> groups = new ArrayList<String>();
//        
//        EqtlSimpleCalc instance = new EqtlSimpleCalc();
//        instance.loadExpression(exprFile);
//        
//        if( _hapMapOnly ){
//            groups.add("CEU");
//            outdir += "_hapmap/";
//            instance.loadHapMapGenotypes(outdir, locus, groups);
//        }
//        else{
//            if( _1000gOnly ){
//                outdir += "_1000g/";
//            }
//            else{
//                outdir += "_1000g_hapmap/";
//                groups.add("CEU");
//            }
//            
//            instance.loadGenotypes(vcfFile, outdir, locus, _1000gOnly, groups);
//        }
//        
//        instance.calcEqtls(outdir);
//    }
    
//    @Test
//    public void performCalc() {
//        ArrayList<String> args = new ArrayList<String>();
//        args.add(Main.OPT_OUTPUT_DIR+"/home/victor/Escritorio/eqtl_mho/");
//        args.add(Main.OPT_LOCUS+"chr5:130600001-136200000");
//        args.add(Main.OPT_VALUE+"hapmap");
//        args.add(Main.OPT_GROUP+"CEU");
//        args.add(Main.OPT_EXPRDATA+"/home/victor/Escritorio/eqtl_mho/LCL_exp_CODE.csv");
//        EqtlSimpleCalc.performCalc(args);
//    }
    
    @Test
    public void calcEqtls2() {
        String testDir = "/home/victor/Escritorio/tests_ngsengine/";
        String tpedFile = testDir+"test_popanalysis.tped";
        String mergedTrkFile = testDir+"test_popanalysis.tracking";
        String outdir = testDir+"eqtl_pca_out/";
        boolean adjustGeno = true;
        
        EqtlSimpleCalc instance = new EqtlSimpleCalc();
        instance.calcEqtls(tpedFile, mergedTrkFile, outdir, adjustGeno);
    }
    
    
}