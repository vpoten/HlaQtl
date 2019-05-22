package org.msgenetics.hlaqtl;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.awsjoblauncher.QueueUtils;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

        
/**
 * Unit test for simple App.
 */
public class MainTest {
    
    public MainTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Main.loadSubjects( new ArrayList<String>() );
        
        String propFile = System.getProperty("user.home")+"/.grails/ngsengine.properties";
        List args = new ArrayList();
        args.add(Main.OPT_PROPS+propFile);
        Main.loadConfig(args);
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
    public void getOptGroups()
    {
        System.out.println("TestGetOptGroups");
        
        ArrayList<String> list = new ArrayList<String>();
        
        List res = (List) Main.getOptGroups(list);
      
        assertNull(res);
        
        list.add("--group=FIN,GBR");
        res = (List) Main.getOptGroups(list);
        
        assertEquals(2, res.size());
    }

    @Test
    public void loadSubjects()
    {
        System.out.println("TestLoadSubjects");
        
        Map subjects = (Map) Main.getSubjects();
        Subject sbj = (Subject) subjects.get("NA11843");
        
        assertEquals(sbj.getId(), "NA11843");
        assertEquals( sbj.getLanes().get(0), "ERR188426" );
        
        sbj = (Subject) subjects.get("NA11881");
        
        assertEquals( sbj.getLanes().get(0), "ERR188272");
        assertEquals( ((List)sbj.getAlleles().get(Main.getDRB5())).get(0), "01");
        assertTrue( ((List)sbj.getAlleles().get(Main.getDRB4())).isEmpty() );
        
        subjects = (Map) Main.getAllSubjects();
        assertEquals(subjects.size(), 91);
    }
    
    
//    @Test
//    public void launchMasterQueue()
//    {
//        System.out.println("TestLaunchMasterQueue");
//        
//        Main.setSqsMaster( QueueUtils.launchMasterQueue(Main.getProperties()) );
//        
//        assertNotNull( Main.getSqsMaster() );
//    }
    
//    @Test
//    public void launchSlaveQueue()
//    {
//        System.out.println("TestLaunchSlaveQueue");
//        
//        Main.setSqsSlave( QueueUtils.launchSlaveQueue(Main.getProperties()) );
//        
//        assertNotNull( Main.getSqsSlave() );
//    }
    
//    @Test
//    public void launchAuditQueues()
//    {
//        System.out.println("TestLaunchAuditQueues");
//        
//        Main.setSqsReceiver( QueueUtils.launchReceiverQueue(Main.getProperties()) );
//        Main.setSqsSender( QueueUtils.launchSenderQueue(Main.getProperties(), 1000*5) );
//        
//        try {
//            Thread.sleep(1000*20);
//        } catch (InterruptedException ex) {
//            Logger.getLogger(MainTest.class.getName()).log(Level.SEVERE, null, ex);
//        }
//      
//        assertNotNull( Main.getSqsSender() );
//        assertTrue( Main.getSqsReceiver().getCount()>1 );
//    }
    
    @Test
    public void getOptSubjects()
    {
        System.out.println("TestGetOptSubjects");
        
        ArrayList list = new ArrayList();
        list.add("--subjs=NA12006");
        
        Map sub = (Map) Main.getOptSubjects( list );
        assertTrue( sub.size()==1 );
        list.clear();
        
        list.add("--subjs=NA12006,NA12045");
        sub = (Map) Main.getOptSubjects( list);
        assertTrue( sub.size()==2);
        assertTrue( sub.containsKey("NA12006") );
    }
    
    
}
