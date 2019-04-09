/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

import tech.tablesaw.api.Table;
import tech.tablesaw.api.ColumnType
import tech.tablesaw.api.DoubleColumn
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.IntColumn
import tech.tablesaw.io.csv.CsvReadOptions


import org.ngsutils.Utils

/**
 *
 * @author victor
 */
class GTExEqtl {
    
    String path
    
    static final regEgenes = /([\w-]+).v7.egenes.txt.gz/
    static final suffEgenes = '.v7.egenes.txt.gz'
    
    // Columns for egenes v7 file:
    static ColumnType[] egenesV7ColumnTypes = [
        ColumnType.STRING, //gene_id:  GENCODE/Ensembl gene ID
        ColumnType.STRING, //gene_name:  GENCODE gene name
        ColumnType.STRING, //gene_chr:  chromosome (gene)
        ColumnType.INTEGER, //gene_start:  gene start position (in base pairs; 1-based coordinates)
        ColumnType.INTEGER, //gene_end:  gene end position (in base pairs; 1-based coordinates)
        ColumnType.STRING, //strand:  genomic strand
        ColumnType.INTEGER, //num_var:  number of variants in cis-window
        ColumnType.DOUBLE, //beta_shape1:  1st shape parameter of the fitted Beta distribution: B(shape1, shape2)
        ColumnType.DOUBLE, //beta_shape2:  2nd shape parameter of the fitted Beta distribution: B(shape1, shape2)
        ColumnType.DOUBLE, //true_df:  Effective degrees of freedom the Beta distribution approximation
        ColumnType.DOUBLE, //pval_true_df
        ColumnType.STRING, //variant_id:  variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37
        ColumnType.INTEGER, //tss_distance:  Distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
        ColumnType.STRING, //chr:  chromosome (variant; same as gene_chr for cis-eQTLs)
        ColumnType.INTEGER, //pos:  position of the first reference base of the variant
        ColumnType.STRING, //ref:  reference sequence of the variant
        ColumnType.STRING, //alt:  alternate sequence of the variant
        ColumnType.INTEGER, //num_alt_per_site:  number of alternative alleles observed at this site
        ColumnType.STRING, //rs_id_dbSNP147_GRCh37p13:  dbSNP142 rsID
        ColumnType.INTEGER, //minor_allele_samples:  number of samples carrying the minor allele
        ColumnType.INTEGER, //minor_allele_count:  total number of minor alleles across individuals
        ColumnType.DOUBLE, //maf:  minor allele frequency observed in the set of donors for a given tissue
        ColumnType.INTEGER, //ref_factor:  '1', when the minor allele is the alt base, '-1' when the minor allele is the reference base
        ColumnType.DOUBLE, //pval_nominal:  nominal p-value associated with the most significant variant for this gene
        ColumnType.DOUBLE, //slope:  regression slope
        ColumnType.DOUBLE, //slope_se:  standard error of the regression slope
        ColumnType.DOUBLE, //pval_perm:  permutation p-value
        ColumnType.DOUBLE, //pval_beta:  beta-approximated permutation p-value
        ColumnType.DOUBLE, //qval:  Storey q-value derived from pval_beta
        ColumnType.DOUBLE, //pval_nominal_threshold:  nominal p-value threshold for calling a variant-gene pair significant for the gene
        ColumnType.DOUBLE, //log2_aFC:  a measure of cis-eQTL effect size, is defined as the log-ratio between the expression of the haplotype carrying the alternative eVariant allele to the one carrying the reference allele.
        ColumnType.DOUBLE, //log2_aFC_lower
        ColumnType.DOUBLE //log2_aFC_upper
    ]
    
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
    Table loadTable(String tissue) {
        def stream = Utils.createInputStream(new File(this.path, "${tissue}${this.suffEgenes}").toString())
        def builder =  CsvReadOptions.builder(stream)
            .tableName(tissue)
            .separator((char)'\t') // table is tab-delimited
            .columnTypes(egenesV7ColumnTypes)

        CsvReadOptions options = builder.build();
        
        Table table = Table.read().usingOptions(options);

        return table
    }
    
    /**
     * Get a table with the best-eqtls from all tissues
     */
    static Table getBestEqtlsAllTissues(String path, Double pvalThr) {
        GTExEqtl instance = new GTExEqtl(path: path)
        
        def filteredTables = instance.getTissues().collect( {tissue -> 
                Table table = instance.loadTable(tissue)
                
                // apply p-value threshold
                def column = (DoubleColumn) table.doubleColumn('qval')
                table = table.where(column.isLessThan(pvalThr))
                
                // add tissue column to the filtered table
                column = StringColumn.create("tissue", table.rowCount())
                (0..table.rowCount()-1).each({column.set(it, tissue)})
                table.addColumns(column)
                return table
            })
        
        def result = filteredTables[0]
        (1..filteredTables.size()-1).each({result.append(filteredTables[it])})
        result.setName("best-eqtls-all-tissues")
        return result
    }
    
    /**
     * Filter eqtls by region 
     */
    static Table filterByRegion(Table srcTable, String chr, int start, int end) {
        def chrCol = (StringColumn) srcTable.stringColumn('chr')
        def snpPosCol = (IntColumn) srcTable.intColumn('pos')
        
        return srcTable.where(chrCol.isEqualTo(chr).and(snpPosCol.isBetweenInclusive(start, end)));
    }
}

