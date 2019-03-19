/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

/**
 * Helper class
 * 
 * @author victor
 */
class CompletionReporter {
	
    double totalJobs = 0.0
    double step = 0.1
    String title = ''
    
    private double doneJobs = 0.0
    private double reported = 0.0
    
    /**
     *
     */
    void update(double value) {
        doneJobs += value
        double rate = doneJobs/totalJobs
        
        if( rate>=(reported+step) ){
            println "${title}: ${rate*100}% completed at ${new Date()}"
            reported = rate
        }
    }
    
    /**
     * thread safe version of update method
     */
    synchronized void updateSafe(double value) {
        update(value)
    }
}

