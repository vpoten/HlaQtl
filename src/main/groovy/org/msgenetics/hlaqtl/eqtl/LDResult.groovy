/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.eqtl

import org.ngsutils.variation.SNPData

/**
 *
 * @author victor
 */
/**
 *
 */
class LDResult {
    SNPData snpBestCis
    SNPData snpTrans
    SNPData snpFeat = null
    Double rSq
}

