/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import tech.tablesaw.api.Table;
import tech.tablesaw.api.ColumnType
import tech.tablesaw.api.DoubleColumn
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.IntColumn

/**
 * Base class for eQTL tables
 * @author victor
 */
class BaseEqtlTable {
    
    static String chrColName = 'chr'
    static String posColName = 'pos'
	
    /**
     * Filter eqtls by region 
     */
    static Table filterByRegion(Table srcTable, String chr, int start, int end) {
        def chrCol = (StringColumn) srcTable.stringColumn(chrColName)
        def snpPosCol = (IntColumn) srcTable.intColumn(posColName)
        
        return srcTable.where(chrCol.isEqualTo(chr).and(snpPosCol.isBetweenInclusive(start, end)));
    }
}

