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
import tech.tablesaw.io.csv.CsvReadOptions

/**
 * Base class for eQTL tables
 * @author victor
 */
abstract class BaseEqtlTable {
    /** path for files */
    String path
    
    /** column name for chr */
    String chrColName = 'chr'
    
    /** column name for snp position */
    String posColName = 'pos'
    
    /** column name for pvalue */
    String pvalueColName = 'pval'
    
    /** column separator */
    Character separator = (char)'\t'
    
    /** get column types */
    abstract ColumnType[] getColumnTypes()
    
    /** get best eqtls types */
    abstract Table getBestEqtls(pvalThr)
    
    /**
     * Load a table from stream
     */
    protected Table _loadTable(stream, name) {
        def builder =  CsvReadOptions.builder(stream)
            .tableName(name)
            .separator(getSeparator()) // table is tab-delimited
            .columnTypes(getColumnTypes())

        CsvReadOptions options = builder.build();
        
        Table table = Table.read().usingOptions(options);

        return table
    }
    
    /**
     * Filter table by p-value threshold
     */ 
    Table filterPvalue(Table table, Double thr) {
        def column = (DoubleColumn) table.doubleColumn(getPvalueColName())
        table = table.where(column.isLessThan(thr))
        return table
    }
	
    /**
     * Filter eqtls by region 
     */
    Table filterByRegion(Table srcTable, String chr, int start, int end) {
        def chrCol = (StringColumn) srcTable.stringColumn(getChrColName())
        def snpPosCol = (IntColumn) srcTable.intColumn(getPosColName())
        
        return srcTable.where(chrCol.isEqualTo(chr).and(snpPosCol.isBetweenInclusive(start, end)));
    }
}

