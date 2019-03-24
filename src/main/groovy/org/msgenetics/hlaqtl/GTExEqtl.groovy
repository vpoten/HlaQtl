/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.CsvReadOptions


import org.ngsutils.Utils

/**
 *
 * @author victor
 */
class GTExEqtl {
    
    static final String path
    
    static final regEgenes = /(\w+).v7.egenes.txt.gz/
    static final suffEgenes = '.v7.egenes.txt.gz'
    
    /**
     * Get tissues availables in path; as extracted from egenes file names
     */
    def getTissues() {
        def egenesFiles = new File(this.path).list({d, f-> f ==~ this.regEgenes} as FilenameFilter).toList()
        return egenesFiles.collect({(it =~ this.regEgenes)[0][1]})
    }
    
    /**
     * Load a table from path by tissue name
     */
    def loadTable(tissue) {
        def stream = Utils.createInputStream(new File(this.path, "${tissue}${this.suffEgenes}").toString())
        def builder =  CsvReadOptions.builder(stream, tissue)
            .separator('\t') // table is tab-delimited

        CsvReadOptions options = builder.build();
        
        Table table = Table.read().usingOptions(options);

        return table
    }
}

