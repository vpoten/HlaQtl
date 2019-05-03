/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table;
import tech.tablesaw.api.StringColumn
import tech.tablesaw.io.csv.CsvReadOptions

import org.ngsutils.Utils

/**
 *
 * @author victor
 */
class IGSRSamples {
    
    /** Internal data table*/
    Table table = null
	
    /**
     * Factory method. Creates an instance loading data from `igsr_samples.tsv` resource file
     */
    static IGSRSamples create() {
        BufferedReader samplesFile = Utils.createReader('igsr_samples.tsv')
        def builder =  CsvReadOptions.builder(samplesFile)
            .tableName('igsr_samples')
            .separator((char)'\t') // table is tab-delimited

        CsvReadOptions options = builder.build();
        
        Table table = Table.read().usingOptions(options);
        
        return new IGSRSamples(table: table)
    }
    
    /**
     * Get subjects codes for the given populations
     */
    List<String> getSubjects(populations) {
        StringColumn column = (StringColumn) table.stringColumn('Population code')
        Table filtered = table.where(column.isIn(populations))
        return filtered.stringColumn('Sample name').asList()
    }
}

