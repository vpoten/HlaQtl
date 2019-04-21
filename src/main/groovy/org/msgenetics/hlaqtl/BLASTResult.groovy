/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl


/**
 *
 */
class BLASTResult {
    List<AlignResult> alignments = []
    String query = null
    String subject = null
    Integer queryLength = null
    Integer subjectLength = null
    
    Integer totalScore() {
        alignments.sum{it.score}
    }
    
    Double averageIdent() {
        Double total = alignments.sum{it.identities}
        return (total==null) ? null : (total/(double)alignments.size())
    }
}

