/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import tech.tablesaw.api.Table
import tech.tablesaw.api.ColumnType

/**
 *
 * @author victor
 */
class SimplifiedEqtl extends BaseEqtlTable {
    
    // Columns for simplified eqtl file:
    static final ColumnType[] simplifiedColumnTypes = [
        ColumnType.STRING, //feat:  feature name or code (gene, transcript, protein, ...)
        ColumnType.STRING, //feat_chr:  chromosome (feature)
        ColumnType.INTEGER, //start:  feature start position (in base pairs; 1-based coordinates)
        ColumnType.INTEGER, //end:  feature end position (in base pairs; 1-based coordinates)
        ColumnType.STRING, //strand:  genomic strand
        ColumnType.STRING, //chr:  chromosome (variant; same as feat_chr for cis-eQTLs)
        ColumnType.INTEGER, //pos:  position of the first reference base of the variant
        ColumnType.STRING, //rs_id: snp rs id
        ColumnType.DOUBLE, //pval:  corrected p-value 
        ColumnType.DOUBLE, //slope:  regression slope
    ]
    
    /**
     * get column types
     */ 
    ColumnType[] getColumnTypes() {
        return simplifiedColumnTypes
    }
    
    /**
     * Load table from path
     */
    Table loadTable(String name) {
        def stream = Utils.createInputStream(new File(this.path).toString())
        return _loadTable(stream, name)
    }
    
    /**
     * get best eqtls
     */
    Table getBestEqtls(pvalThr) {
        Table table = loadTable('best_eqtls')
        return filterPvalue(table, pvalThr)
    }
}
