/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl.hla

import org.biojava.nbio.core.sequence.AccessionID
import org.biojava.nbio.core.sequence.DNASequence
import org.biojava.nbio.core.sequence.Strand
import org.biojava.nbio.core.sequence.compound.DNACompoundSet
import org.biojava.nbio.alignment.Alignments
import org.biojava.nbio.alignment.SubstitutionMatrixHelper
import org.biojava.nbio.alignment.SimpleSubstitutionMatrix
import org.biojava.nbio.alignment.template.SubstitutionMatrix
import org.biojava.nbio.alignment.SimpleGapPenalty
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType
import org.biojava.nbio.core.sequence.template.Sequence
import org.biojava.nbio.core.sequence.compound.NucleotideCompound
import org.biojava.nbio.core.sequence.compound.DNACompoundSet
import org.biojava.bio.program.gff3.GFF3Record
import org.biojava.bio.program.gff.GFFRecord
import org.ngsutils.annotation.*
import org.ngsutils.annotation.biojava.SimpleGFFRecordLight
import org.ngsutils.annotation.biojava.GFF3DocumentHandlerConvert as GFF3D

/**
 *
 * @author victor
 */
class HLAAllele {
    
    String name
    String hlaId
    DNASequence cds
    DNASequence geneSeq = null
    boolean isPartial = false
    HLAAllele parent // most similar allele with complete gene sequence
    Gene gene
    
    /**
     * 
     */ 
    int getAnnotationLength() {
        Isoform iso = gene.getIsoforms().get(name)
        iso.exons.sum{ it.end-it.start+1 }
    }
    
    /**
     *
     */
    def writeAnnotation(int offset, String seqId, Writer writer) {
        Isoform rna = gene.isoforms.get(name)
        def idSuffix = '_'+seqId
        
        def genemap = [:]
        //write gene record
        genemap.seqid = seqId
        genemap.source = 'hlaqtl'
        genemap.type = 'gene'
        genemap.start = (rna.exons[0].start as Integer) + offset
        genemap.end = (rna.exons[rna.exons.size()-1].end as Integer) + offset
        genemap.score = '.'
        genemap.strand = rna.exons[0].strand.token
        genemap.phase = '.'
        genemap.attributes = ["${GFF3D.TERM_ID}":gene.geneId+idSuffix]
        GFF3D.writeGff3Record(writer, genemap)
        writer.write("\n")

           
        //write mrna
        def mrna = [:]

        mrna.seqid = seqId
        mrna.source = 'hlaqtl'
        mrna.type = 'mRNA'
        mrna.start = (rna.exons[0].start as Integer) + offset
        mrna.end = (rna.exons[rna.exons.size()-1].end as Integer) + offset
        mrna.score = '.'
        mrna.strand = rna.exons[0].strand.token
        mrna.phase = '.'
        mrna.attributes = ["${GFF3D.TERM_ID}":rna.exons[0].transcriptId+idSuffix, "${GFF3D.TERM_PARENT}":rna.exons[0].geneId+idSuffix]
        GFF3D.writeGff3Record(writer, mrna)

        writer.write("\n")

        //write CDSs (exons)
        rna.exons.each{ rec->
            def cds = [:]

            cds.seqid = seqId
            cds.source = 'hlaqtl'
            cds.score = (rec.score==GFFRecord.NO_SCORE) ? '.' : rec.score
            cds.type = 'CDS'
            cds.start = rec.start + offset
            cds.end = rec.end + offset
            cds.phase = (rec.frame==GFFRecord.NO_FRAME) ? '0' : rec.frame
            cds.strand = rec.strand.token
            cds.attributes = ["${GFF3D.TERM_PARENT}":rec.transcriptId+idSuffix]
            GFF3D.writeGff3Record(writer, cds)
        }

        writer.write("\n\n")
    }
    
    /**
     * build gene sequence and simultaneously rebuilds the annotation
     */
    def buildGeneSeq() {
        if( geneSeq!=null ){
            System.err.println("HLAAllele: gene sequence/annotation already rebuilded.")
            return
        }
        
        Isoform iso = gene.getIsoforms().get(name)
        Isoform isoParent = gene.getIsoforms().get(parent.name)
        
        StringBuilder newSeq = new StringBuilder()
        int seqPos = 0 //zero based
        String geneStr = parent.geneSeq.getSequenceAsString()
        def newExonAnnotation = []
        
        // use annotation and geneSeq of parent to build the sequence
        isoParent.getExons().each{ exon->
            SimpleGFFRecordLight modExon = iso.getExons().find{ it.exonNumber==exon.exonNumber }
            
            SimpleGFFRecordLight newExon = new SimpleGFFRecordLight(modExon ?: exon)
            newExon.isoform = iso
            newExon.transcriptId = iso.transcriptId
            
            String exonStr = parent.geneSeq.getSequenceAsString(exon.start, exon.end, Strand.POSITIVE)
            int actualStart = geneStr.indexOf( exonStr )
            
            //add intron/UTR to sequence
            newSeq << geneStr.substring(seqPos, actualStart)
            
            newExon.start = newSeq.length()+1
            
            if( modExon ){
                try{
                newSeq << this.cds.getSequenceAsString(modExon.start, modExon.end, Strand.POSITIVE)
                } catch(e){
                    println "Error bulding seq for gene: ${name}, parent ${isoParent.transcriptId}"
                    throw e
                }
            }
            else{
                newSeq << exonStr
            }
            
            newExon.end = newSeq.length()
            
            seqPos = actualStart + exonStr.length()
            newExonAnnotation << newExon
        }
        
        //add last UTR
        if( seqPos<geneStr.length() ) {
            newSeq << geneStr.substring(seqPos)
        }
        
        iso.exons = newExonAnnotation
        geneSeq = new DNASequence(newSeq.toString())
        geneSeq.accession = new AccessionID("HLA:${hlaId} ${name} ${newSeq.length()} bp")
    }
    
    /**
     * rebuild annotation for alleles with complete sequence in _gen.fasta
     */
    def rebuildAnnotation() {
        // Nothing to do : seems to be corrected in current release of HLA db.
    }
    
    /**
     *
     */
    def buildGeneSeqFromIntrons(List introns) {
        StringBuilder newSeq = new StringBuilder()
        def newExonAnnotation = []
        
        Isoform isoParent = gene.getIsoforms().get(parent ? parent.name : name)
        Isoform iso = gene.getIsoforms().get(name)
        
        isoParent.getExons().eachWithIndex{ exon, i->
            SimpleGFFRecordLight modExon = iso.getExons().find{ it.exonNumber==exon.exonNumber }
            
            SimpleGFFRecordLight newExon = new SimpleGFFRecordLight(modExon ?: exon)
            newExon.isoform = iso
            newExon.transcriptId = iso.transcriptId
            
            newExon.start = newSeq.length()+1
            
            if( modExon ){
                try{
                newSeq << this.cds.getSequenceAsString(modExon.start, modExon.end, Strand.POSITIVE)
                } catch(e){
                    println "Error bulding seq for gene: ${name}, parent ${isoParent.transcriptId}"
                    throw e
                }
                
            }
            else{
                newSeq << parent.geneSeq.getSequenceAsString(exon.start, exon.end, Strand.POSITIVE)
            }
            
            newExon.end = newSeq.length()
            
            if( i<introns.size() ) {
                newSeq << introns[i].getSequenceAsString()
            }
            
            newExonAnnotation << newExon
        }
        
        iso.exons = newExonAnnotation
        geneSeq = new DNASequence(newSeq.toString())
        geneSeq.accession = new AccessionID("HLA:${hlaId} ${name} ${newSeq.length()} bp")
    }
    
    /**
     * 
     */ 
    DNASequence getCompleteCDS() {
        if( !isPartial ){
            return this.cds
        }
        
        StringBuilder strb = new StringBuilder()
        Isoform iso = gene.getIsoforms().get(name)
        
        iso.exons.each{ exon->
            strb << this.geneSeq.getSequenceAsString(exon.start, exon.end, Strand.POSITIVE)
        }
        
        return new DNASequence(strb.toString(), DNACompoundSet.getDNACompoundSet())
    }
    
    /**
     *
     */
    def alignment(Sequence target, boolean isLocal, boolean useGeneSeq=false) {
        SubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix(
                DNACompoundSet.getDNACompoundSet(), 2 as Short, -6 as Short );
 
        SimpleGapPenalty gapP = new SimpleGapPenalty(5 as Short, 2 as Short);

        PairwiseSequenceAligner<Sequence, NucleotideCompound> psa =
                        Alignments.getPairwiseAligner((useGeneSeq ? this.geneSeq : this.completeCDS), 
                            target, (isLocal ? PairwiseSequenceAlignerType.LOCAL : PairwiseSequenceAlignerType.GLOBAL), 
                            gapP, matrix);
                                    
        return psa
    }
    
}

