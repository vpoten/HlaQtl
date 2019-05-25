/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import org.awsjoblauncher.postprocess.CommConstants;
import org.awsjoblauncher.postprocess.JobDataBean;
import org.awsjoblauncher.postprocess.PostBase;
import org.awsjoblauncher.storage.StorageManager;
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
public class JobGeneratorTest {
    
    static Subject subject = new Subject();
    
    public JobGeneratorTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Properties props=org.awsjoblauncher.Main.loadProperties("ngsengine.properties");
        PostBase.setProps(props);
        
        subject.setId("NA11881");
        subject.getLanes().add("1982_8");
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
    public void fastQC() {
        System.out.println("TestFastQC");
        
        String inputUrl = "s3://raw_data/ceu/";
        
        JobDataBean data = JobGenerator.fastQC(subject, inputUrl);
        
        assertNotNull(data);
        
        String outdir = StorageManager.OUTPUT_PRE+CommConstants.COMMNAME_FASTQC+"_"+subject.getId();
        String inputs = inputUrl+subject.getLanes().get(0)+"_1.fastq.gz ";
        inputs += inputUrl+subject.getLanes().get(0)+"_2.fastq.gz ";
        
        assertTrue( data.getCommand().endsWith(
                CommConstants.FASTQC_FOLDER+"fastqc -q -t 2 -o "+outdir+" "+inputs
                ) );
    }
    
    @Test
    public void fastxToolkit() {
        System.out.println("TestFastxToolkit");
        
        String inputUrl = "s3://raw_data/ceu/";
        
        //{minLen, maxLen, trimFirst, qualFilt, percQual}
        HashMap optMap = new HashMap();
        optMap.put("minLen", "10");
        optMap.put("maxLen", "200");
        optMap.put("trimFirst", "12");
        optMap.put("qualFilt", "20");
        optMap.put("percQual", "80");
        
        JobDataBean data = JobGenerator.fastxToolkit(subject, inputUrl, optMap);
        
        assertNotNull(data);
    }
    
    @Test
    public void tophat() {
        System.out.println("TestTophat");
        
        String input = "output://fastx_toolkit_NA11881/";
        JobDataBean data = JobGenerator.tophat(subject, input, "/index/hg19");
        assertNotNull(data);
    }
    
    @Test
    public void cufflinks() {
        System.out.println("TestCufflinks");
        
        String input = "output://tophat_NA11881/";
        JobDataBean data = JobGenerator.cufflinks(subject, input, null, false);
        assertNotNull(data);
        
        data = JobGenerator.cufflinks(subject, input, "annotations://9606/refGene.gtf", false);
        assertNotNull(data);
        assertTrue( data.getUrl().contains("_ref_NA11881") );
        
        data = JobGenerator.cufflinks(subject, input, "annotations://9606/refGene.gtf", true);
        assertNotNull(data);
        assertTrue( data.getUrl().contains("_refOnly_NA11881") );
        
        data = JobGenerator.cufflinks(subject, input, "annotations://9606/ensGene.gtf", true);
        assertTrue( data.getUrl().contains("_ensOnly_NA11881") );
    }
    
    @Test
    public void cuffcompare() {
        System.out.println("TestCuffcompare");
        
        String input = "cufflinks://tophat_NA11881/";
        ArrayList list = new ArrayList();
        list.add(input);
        JobDataBean data = JobGenerator.cuffcompare(subject, list, "annotations://9606/ref.gtf");
        assertNotNull(data);
    }
    
}
