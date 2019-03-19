/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.msgenetics.hlaqtl

/**
 *
 */
class AlignResult {
    Integer score = null
    Integer scoreBits = null
    Double expect = null
    Double identities
    Double gaps
    String strand
}

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

/**
 *
 * @author victor
 */
class SimpleBLAST {
    
    static final String NO_HITS = '***** No hits found *****'
    static final String SCORE = 'Score'
    static final String IDENT = 'Identities'
    static final String STRAND = 'Strand='
    static final String QUERY = 'Query='
    static final String SUBJECT = 'Subject='
    static final String LENGTH = 'Length='
    
    String binPath //ncbi BLAST+ binaries
    
    /**
     *
     */
    public SimpleBLAST(path) {
        binPath = path
        
        if( !binPath.endsWith('/') ){
            binPath += '/'
        }
    }
    
    /**
     *
     */
    public BLASTResult runBlastn(String query, String subject) {
        def proc = "${binPath}blastn -query ${query} -subject ${subject}".execute()
       
        def getRate = { str->
            def toks = str.split("/")
            return ((toks[0] as Double) / (toks[1] as Double))
        }
        
        def removeBrackets = { str->
            return str.substring(str.indexOf('(')+1, str.indexOf(')'))
        }
        
        // parse blasn result
        def reader =  new BufferedReader( new InputStreamReader(proc.getInputStream()) )
        
        String line = null
        BLASTResult res = new BLASTResult()
        AlignResult current = null
        
        while( (line = reader.readLine())!=null ){   
            line = line.trim()
            
            if( line.isEmpty() )
                continue;
            
            if( line==NO_HITS ){
                return new BLASTResult()
            }
            else if( line.startsWith(QUERY) ){
                res.query = line.substring(QUERY.length()).trim()
            }
            else if( line.startsWith(SUBJECT) ){
                res.subject = line.substring(SUBJECT.length()).trim()
            }
            else if( line.startsWith(LENGTH) ){
                Integer length = line.substring(LENGTH.length()).trim() as Integer
                
                if( res.queryLength==null ){
                    res.queryLength = length
                }
                else{
                    res.subjectLength = length
                }
            }
            else if( line.startsWith(SCORE) ){
                current = new AlignResult()
                def toks = (line.split("\\s") as List).findAll{it}
                current.scoreBits = (toks[2] as Double) as Integer//score bits
                current.score = (removeBrackets(toks[4]) as Double) as Integer//score
                current.expect = toks[7] as Double //expect
            }
            else if( line.startsWith(IDENT) ){
                def toks = (line.split("\\s") as List).findAll{it}
                current.identities = getRate(toks[2]) //identities
                current.gaps = getRate(toks[6]) //gaps
            }
            else if( line.startsWith(STRAND) ){
                current.strand = line.substring(STRAND.length()).trim()
                res.alignments << current
            }
        }
        
        reader.close()
        
        try {
            //wait for process termination
            proc.waitFor()
        } catch (InterruptedException ex) {
        }
        
        if( proc.exitValue()!=0 ){
            // error
            return null
        }
                
        return res
    }
    
}

