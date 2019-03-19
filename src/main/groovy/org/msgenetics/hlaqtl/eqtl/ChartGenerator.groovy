/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.ngsutils.eqtl.CorrelationCalc

import org.jfree.chart.*
import org.jfree.data.category.DefaultCategoryDataset
import org.jfree.data.xy.DefaultXYDataset
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset
import org.jfree.chart.plot.PlotOrientation
import org.jfree.chart.renderer.xy.XYDotRenderer
import org.jfree.chart.plot.ValueMarker
import org.jfree.chart.axis.NumberAxis
import org.jfree.ui.Layer
import org.jfree.chart.renderer.category.MinMaxCategoryRenderer
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer
import java.awt.Color

/**
 *
 * @author victor
 */
class ChartGenerator {
    CorrelationCalc correlation
    String outDir
    def snps //map of {snpid, SNPData object}
    
    //constants
    static final int CHART_WIDTH = 640
    static final int CHART_HEIGHT = 480
    static final String CHART_PREF = 'chart_'
    static final double COORD_CONV = 1e-6
    static final double EXPR_MIN_VAL = 3.0
    
    /**
     * generate associated charts
     */
    def generate(double pvalueThr, double corrThr, int total){
        def selIsoforms = [] as TreeSet
        
        int total1 = 0
        
        //select snps and isoforms
        correlation.snps.each{ snp->
            correlation.isoforms.each{ iso->
                Double pval = correlation.getCorrPValues(snp,iso)
                Double corr = correlation.getCorrValues(snp,iso)
                
                if(  pval!=null && corr!=null && pval<pvalueThr ){
                    corr = Math.abs(corr)
                    
                    if( corr > corrThr ){
                        if( genotypeExpres(snp, iso) ){
                            println "Generate genotype express. chart for ${snp} - ${iso} [pval=${pval}, corr=${corr}]"
                            total1++
                            selIsoforms << iso
                        }
                        else{
                           println "Discarded genotype express. chart for ${snp} - ${iso} [pval=${pval}, corr=${corr}]" 
                        }
                    }
                }
            }
        }
        
        selIsoforms.each{
            println "Generate snp assoc. region chart for ${it}"
            snpAssocRegion(it)
        }
        
        println "${total1} genotype express. charts generated"
        println "${selIsoforms.size()} snp assoc. region charts generated"
    }
    
    /**
     *
     */
    private double[][] toDoubleArray(list){
        def array = new double [list.size()][list[0].size()]
        list.eachWithIndex{ l2, i->
            l2.eachWithIndex{ ele, j->
                array[i][j] = ele
            }
        }
        return array
    }
    
    /**
     * generates scatter plot: X=snp coordinates, Y=correlation (abs)
     */
    private def snpAssocRegion(isoform){
        def dataset = new DefaultXYDataset()
        def datarows = [[],[]]
        
        snps.each{ snpId, snpData->
            Double val = correlation.getCorrValues(snpId, isoform)
            
            if( val!=null ){
                datarows[0] << snpData.position*COORD_CONV
                datarows[1] << Math.abs(val)
            }
        }
        
        // add value to Dataset
        dataset.addSeries( isoform, toDoubleArray(datarows))
        
        // render chart
        def imgFile = "${CHART_PREF}snpreg_${isoform}.png"

        exportScatterPlot(dataset, "${isoform} SNPs association", 
            'Genomic coordinate', 'Correlation', outDir+imgFile )
    }
    
    /**
     * generates box and whisker plot: X=genotypes, Y=expression
     */
    private def genotypeExpres(snp, isoform){
        def dataset = new DefaultBoxAndWhiskerCategoryDataset()
        def snpData = snps[snp]
        
        def exprVals = [0:[],1:[],2:[]]
        
        correlation.subjects.each{ sbjId->
            Double genot = correlation.genotypes[sbjId][snp]
            Double expr = correlation.expression[sbjId][isoform]
            
            if(genot!=null && expr!=null){
                exprVals[(int)genot] << expr
            }
        }
        
        // discard chart if max expr. value < MIN_VAL
        def max = [0,1,2].collect{exprVals[it].max()}.max()
        if( max < EXPR_MIN_VAL ){
            return false
        }
        
        [0,1,2].each{
            //add list to dataset: category=genotype, series=expres.
            dataset.add( exprVals[it], 'Expression', snpData.decode(it) )
        }
        
        // render chart
        def imgFile = "${CHART_PREF}genoexpr_${snp}_${isoform}.png"

        exportBoxChart(dataset, "${snp} ${isoform} genotype expression", 
            'Genotype', 'Expression', outDir+imgFile )
        
        return true
    }
    
    /**
     * writes chart to disk (PNG)
     */
    private static exportScatterPlot(xydataset, title, xAxisLabel, yAxisLabel, String outName){
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, xAxisLabel, yAxisLabel, 
            xydataset, PlotOrientation.VERTICAL, true, true, false)
        
        def xyplot = jfreechart.getPlot()
        xyplot.setDomainCrosshairVisible(true)
        xyplot.setDomainCrosshairLockedOnData(true)
        xyplot.setRangeCrosshairVisible(true)
        xyplot.setRangeCrosshairLockedOnData(true)
        xyplot.setDomainZeroBaselineVisible(true)
        xyplot.setRangeZeroBaselineVisible(true)
        XYDotRenderer xydotrenderer = new XYDotRenderer()
        xydotrenderer.setDotWidth(2)
        xydotrenderer.setDotHeight(2)
        xyplot.setRenderer(xydotrenderer)
        def numberaxis = xyplot.getDomainAxis()
        numberaxis.setAutoRangeIncludesZero(false)
        
        ChartUtilities.saveChartAsPNG( new File(outName), jfreechart, CHART_WIDTH, CHART_HEIGHT)
    }
    
    /**
     * writes chart to disk (PNG)
     */
    private static exportBoxChart(dataset, title, categoryAxisLabel, valueAxisLabel, String outName){
        JFreeChart jfreechart =  ChartFactory.createBoxAndWhiskerChart( title,
              categoryAxisLabel, valueAxisLabel, dataset, false)
         
        def categoryplot = jfreechart.getPlot()
        categoryplot.setDomainGridlinesVisible(true)
        categoryplot.setRangePannable(true)
        NumberAxis numberaxis = (NumberAxis)categoryplot.getRangeAxis()
        numberaxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits())
        
        def renderer = new BoxAndWhiskerRenderer()
        renderer.setMeanVisible(false)
        renderer.setMedianVisible(true)
        categoryplot.setRenderer(renderer)
        
        ChartUtilities.saveChartAsPNG( new File(outName), jfreechart, CHART_WIDTH, CHART_HEIGHT)
    }
    
}

