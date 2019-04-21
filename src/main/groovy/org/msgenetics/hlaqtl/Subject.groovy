/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

/**
 *
 */
class Subject{
    String id
    
    // map of keys = {DRA*,DRB1*,..} and values = list of haplotypes codes (01:10, 01, ...)
    Map alleles = [:]
    
    // string coded alleles
    List haplotypes = ['','']
    
    List lanes = [] // list of lane codes
    List fastqReads = [] // list of FastqReads
    
    /**
     * encodes haplotypes using alleles map info.
     */
    def encodeHaplotypes(){
        Main.hlaGenes.each{ gene->
            def list = this.alleles[gene]
            (0..1).each{ haplotypes[it] += ',' }

            if( list.size()==1 ){
                (0..1).each{ haplotypes[it] += (gene+list[0]) }
            }
            else if( list.size()==2 ){
                (0..1).each{ haplotypes[it] += (gene+list[it]) }
            }
            else{
                (0..1).each{ haplotypes[it] += gene }
            }
        }
        
        (0..1).each{ haplotypes[it] = haplotypes[it].substring(1) }
        haplotypes = haplotypes.sort()
    }
    
    /**
     *
     */
    String getGenotype() {
        return "${haplotypes[0]} ${haplotypes[1]}"
    }
    
}

