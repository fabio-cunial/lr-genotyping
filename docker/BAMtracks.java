import java.util.Arrays;
import java.util.HashMap;
import java.io.*;


/**
 * 
 */
public class BAMtracks {
    
    private static final String SAM_SEPARATOR = "\t";
    
	/**
	 * Indels shorter than this are ignored
	 */
	private static final int SV_LENGTH = 40;
    
    /**
     * Clips shorter than this do not create events
     */
    private static final int MIN_CLIP_LENGTH = SV_LENGTH;
    
	/**
	 * SAM records with QUAL field smaller than this are discarded. In practice 
	 * the distribution of QUAL values peaks at 0 and 60, and it has very small 
	 * mass between these values.
	 */
	private static final int QUAL_MIN = 1;
	
	/**
	 * Secondary alignments are usually repeat-induced and are likely discarded
	 * by SV callers.
	 */
	private static final boolean DISCARD_SECONDARY = true;
	
	/**
	 * Supplementary alignments should be components of alignment chains and are
	 * not likely to be discarded by SV callers.
	 */
	private static final boolean DISCARD_SUPPLEMENTARY = false;
	
	/**
	 * Sliding window on the reference, over which the tracks are computed.
	 */
	private static int WINDOW_LENGTH, WINDOW_STEP;
    
    /**
     * For comparing INS sequences
     */
    private static int KMER_LENGTH;
    private static final boolean CANONIZE_KMERS = false;
    private static int MAX_KMER_VECTORS = 1000;  // Arbitrary
    
    /**
     * Event types
     */
    private static final int DEL_START = 0;
    private static final int DEL_END = 1;
    private static final int INS_POS = 2;
    private static final int CLIP_LEFT = 3;  // The alignment is on the left
    private static final int CLIP_RIGHT = 4;  // The alignment is on the right
	
	/**
	 * All the BAM records that intersect the current window
	 */
	private static Alignment[] alignments;
	private static int lastAlignment;
    
    /**
     * All events in the current window
     */
    private static int[] del, ins, clip;
    private static int del_last, ins_last, clip_last;
    private static Kmer[][] kmerVectors;
    private static int[] kmerVectors_alignmentID;
    private static int kmerVectors_last;
	
	/**
	 * Temporary space
	 */
	private static int[] tmpArray = new int[3];
    private static StringBuilder sb;
    private static HashMap<Kmer,Kmer> map;
    private static Kmer[] kmerPool, kmerVector;
    private static int[] coverageHistogram;
	
	
	/**
	 *
	 */
	public static void main(String[] args) throws IOException {
		int i, p, q;
		int contig, flags, quality, position, currentContig, currentStart, nAlignments;
        double rate;
		String str, cigar, seq;
		BufferedReader br;
		BufferedWriter bw;
		
		// Parsing the input
		final String SAM_FILE = args[0];  // Sorted
		final String OUTPUT_FILE = args[1];
        WINDOW_LENGTH=Integer.parseInt(args[2]);
        WINDOW_STEP=Integer.parseInt(args[3]);
        KMER_LENGTH=Integer.parseInt(args[4]);
		
		// Allocating memory
		alignments = new Alignment[50];  // Arbitrary
		for (i=0; i<alignments.length; i++) alignments[i] = new Alignment();
        del = new int[400];  // Arbitrary
        ins = new int[300];  // Arbitrary
        clip = new int[300];  // Arbitrary
        kmerVectors = new Kmer[100][0];  // Arbitrary
        kmerVectors_alignmentID = new int[100];  // Arbitrary
        sb = new StringBuilder();
        map = new HashMap<Kmer,Kmer>();
        kmerPool = new Kmer[100];  // Arbitrary
        for (i=0; i<kmerPool.length; i++) kmerPool[i] = new Kmer();
        kmerVector = new Kmer[100];  // Arbitrary
        coverageHistogram = new int[WINDOW_LENGTH];
		
		// Scanning the reference
		br = new BufferedReader(new FileReader(SAM_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		nAlignments=0; currentContig=1; currentStart=1; lastAlignment=-1; 
		str=br.readLine();
		while (str!=null) {
			nAlignments++;
			p=str.indexOf(SAM_SEPARATOR); q=str.indexOf(SAM_SEPARATOR,p+1);
			flags=Integer.parseInt(str.substring(p+1,q));
			if ((DISCARD_SECONDARY&&((flags&0x100)!=0)) || (DISCARD_SUPPLEMENTARY&&((flags&0x800)!=0)) || (flags&0x200)!=0 || (flags&0x400)!=0) {
				str=br.readLine();
				continue;
			}
			p=q; q=str.indexOf(SAM_SEPARATOR,p+1);
			contig=string2contig(str.substring(p+1,q));
			p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
			position=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
			quality=Integer.parseInt(str.substring(p,q));
			if (quality<QUAL_MIN) {
				str=br.readLine();
				continue;
			}
			if (contig!=currentContig) {
                getTracks(currentContig,currentStart,bw,CANONIZE_KMERS);
				lastAlignment=-1; currentContig=contig; currentStart=1;
			}
			else {
                while (position>currentStart+WINDOW_LENGTH-1) {
    				getTracks(currentContig,currentStart,bw,CANONIZE_KMERS);
    				currentStart+=WINDOW_STEP;
    				filterAlignments(currentStart);
                }
			}
            p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
            cigar=str.substring(p,q);
            p=str.indexOf(SAM_SEPARATOR,q+1);
            p=str.indexOf(SAM_SEPARATOR,p+1);
            p=str.indexOf(SAM_SEPARATOR,p+1);
            q=str.indexOf(SAM_SEPARATOR,p+1);
            seq=str.substring(p+1,q);
			appendAlignment(position,cigar,seq);
			if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments");
			str=br.readLine();
		}
        getTracks(currentContig,currentStart,bw,CANONIZE_KMERS);
		br.close(); bw.close();
	}

	
	private static final void appendAlignment(int position, String cigar, String seq) {
		int i;
		final int alignmentsLength = alignments.length;
		
		lastAlignment++;
		if (lastAlignment==alignmentsLength) {
			Alignment[] newArray = new Alignment[alignmentsLength<<1];
			System.arraycopy(alignments,0,newArray,0,alignmentsLength);
			for (i=alignmentsLength; i<newArray.length; i++) newArray[i] = new Alignment();
			alignments=newArray;
		}
		alignments[lastAlignment].set(position,cigar,seq);
	}
	
	
	/**
	 * Removes from $alignments$ every alignment whose projection on the
	 * reference does not intersect $[from..from+WINDOW_LENGTH-1]$ (1-based).
	 */
	private static final void filterAlignments(int from) {
		int i, j;
		final int to = from+WINDOW_LENGTH-1;
		Alignment tmpAlignment;
		
		j=-1;
		for (i=0; i<=lastAlignment; i++) {
			if (Math.min(to,alignments[i].endA)>=Math.max(from,alignments[i].startA)) {
				j++;
				tmpAlignment=alignments[i];
				alignments[i]=alignments[j];
				alignments[j]=tmpAlignment;
			}
		}
		lastAlignment=j;
	}
    
    
    /**
     * Given window $currentContig[currentStart..currentStart+WINDOW_LENGTH-1]$,
     * for each event type, the procedure outputs the max distance between an 
     * event in the window, and its nearest neighbor on a different alignment in
     * the window. The procedure outputs 0 iff there is no event of a given type 
     * in the window; it outputs -1 iff there is an event with no neighbor on a
     * different alignment.
     */
    private static final void getTracks(int currentContig, int currentStart, BufferedWriter bw, boolean canonizeKmers) throws IOException {
        boolean isolated;
        int i, j;
        int dx, dy;
        double sum, distance, minDistance, maxDel, maxIns, maxClip, maxKmer, avgCoverage;
        
        if (lastAlignment==-1) {
            bw.write(currentContig+","+currentStart+",0,0,0,0,0\n");
            return;
        }
        
        // Collecting events (in no particular order).
        del_last=-1; ins_last=-1; clip_last=-1; kmerVectors_last=-1;
        Arrays.fill(coverageHistogram,0);
        for (i=0; i<=lastAlignment; i++) alignments[i].addEvents(currentStart,i,canonizeKmers,coverageHistogram);
        avgCoverage=0.0;
        for (i=0; i<WINDOW_LENGTH; i++) avgCoverage+=coverageHistogram[i];
        avgCoverage/=WINDOW_LENGTH;
        
        // DEL
        maxDel=0.0;
        for (i=0; i<=del_last; i+=4) {
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=del_last; j+=4) {
                if (j==i || del[j+3]==del[i+3] || del[j]!=del[i]) continue;
                dx=del[j+1]-del[i+1]; dy=del[j+2]-del[i+2];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  // Isolated event
                maxDel=-1;
                break;
            }
            else if (minDistance>maxDel) maxDel=minDistance;
        }
        // INS
        maxIns=0.0;
        for (i=0; i<=ins_last; i+=3) {
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=ins_last; j+=3) {
                if (j==i || ins[j+2]==ins[i+2]) continue;
                dx=ins[j]-ins[i]; dy=ins[j+1]-ins[i+1];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  // Isolated event
                maxIns=-1;
                break;
            }
            else if (minDistance>maxIns) maxIns=minDistance;
        }
        // CLIP
        maxClip=0.0;
        for (i=0; i<=clip_last; i+=3) {
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=clip_last; j+=3) {
                if (j==i || clip[j+2]==clip[i+2] || clip[j]!=clip[i]) continue;
                distance=clip[j+1]-clip[i+1];
                if (distance<0) distance=-distance;
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  // Isolated event
                maxIns=-1;
                break;
            }
            else if (minDistance>maxClip) maxClip=minDistance;
        }
        // KMERS
        if (kmerVectors_last+1>MAX_KMER_VECTORS) maxKmer=-2; 
        else {
            for (i=0; i<=kmerVectors_last; i++) {
                sum=0.0;
                for (j=0; j<kmerVectors[i].length; j++) sum+=kmerVectors[i][j].count;
                for (j=0; j<kmerVectors[i].length; j++) kmerVectors[i][j].count/=sum;
            }
            maxKmer=0.0;
            for (i=0; i<=kmerVectors_last; i++) {
                minDistance=Integer.MAX_VALUE;
                for (j=0; j<=kmerVectors_last; j++) {
                    if (j==i || kmerVectors_alignmentID[j]==kmerVectors_alignmentID[i]) continue;
                    distance=kmerDistance(i,j);
                    if (distance<minDistance) minDistance=distance;
                }
                if (minDistance==Integer.MAX_VALUE) {  // Isolated event
                    maxKmer=-1;
                    break;
                }
                else if (minDistance>maxKmer) maxKmer=minDistance;
            }
        }
        
        bw.write(currentContig+","+currentStart+","+maxDel+","+maxIns+","+maxClip+","+maxKmer+","+avgCoverage+"\n");
    }
    
    
    
    
    // ------------------------- CIGAR PROCEDURES ------------------------------
	
	private static class Alignment {
		/**
		 * CIGAR string, and projection of the alignment on it (0-based).
		 */
        public String cigar;
		public int startCigar, endCigar;
		
		/**
		 * Projection of the alignment on the contig (A), 1-based.
		 */
		public int startA, endA;
		
		/**
		 * Substring (B) of the read, possibly soft-clipped, and projection of
         * the alignment on it (0-based).
		 */
        public String seq;
		public int startB, endB;
        
        /**
         * Sum of hard and soft clips
         */
        public int leftClip, rightClip;
		
		
		public final void set(int position, String cigar, String seq) {
			this.cigar=cigar; this.seq=seq.toLowerCase();
            startA=position;
            cigar_prefixClip(cigar,tmpArray);
			startB=tmpArray[0]; startCigar=tmpArray[1]; leftClip=tmpArray[2];
			cigar_length(cigar,startCigar,tmpArray);
			endCigar=tmpArray[0]; endA=startA+tmpArray[1]-1; endB=startB+tmpArray[2]-1;
            cigar_suffixClip(cigar,tmpArray);
            rightClip=tmpArray[0];
		}
		
        
		/**
         * Adds to the global lists all alignment events that fall in the window
         * $[referenceFrom .. referenceFrom+WINDOW_LENGTH-1]$. The procedure 
         * ignores mismatches, and indels shorter than $SV_LENGTH$.
		 */
		public final void addEvents(int referenceFrom, int alignmentID, boolean canonizeKmers, int[] coverageHistogram) {
			char c;
			int i, j, p;
			int length, first, last, referencePosition, seqPosition;
			final int referenceTo = referenceFrom+WINDOW_LENGTH-1;
            
            // CLIP events
            if (startA>=referenceFrom && startA<=referenceTo && leftClip>=MIN_CLIP_LENGTH) addEvent(CLIP_RIGHT,startA,0,alignmentID);
            if (endA>=referenceFrom && endA<=referenceTo && rightClip>=MIN_CLIP_LENGTH) addEvent(CLIP_LEFT,endA,0,alignmentID);
            
            // INS and DEL events
			referencePosition=startA-1; seqPosition=startB-1; p=startCigar;
			for (i=startCigar; i<=endCigar; i++) {
				c=cigar.charAt(i);
				if (c=='M' || c=='X' || c=='=') {
					length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
					for (j=0; j<length; j++) {
						referencePosition++; seqPosition++;
						if (referencePosition>referenceTo) break;
						else if (referencePosition<referenceFrom) continue;
                        coverageHistogram[referencePosition-referenceFrom]++;
					}
				}
				else if (c=='D' || c=='N') {
					length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                    first=referencePosition+1; last=referencePosition+length;
                    if (length>=SV_LENGTH) {
                        if (first>=referenceFrom && first<=referenceTo) addEvent(DEL_START,first,length,alignmentID);
                        if (last>=referenceFrom && last<=referenceTo) addEvent(DEL_END,last,length,alignmentID);
                    }
					for (j=0; j<length; j++) {
						referencePosition++;
						if (referencePosition>referenceTo) break;
						else if (referencePosition<referenceFrom) continue;
					}
				}
                else if (c=='I') {
                    length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                    if (length<SV_LENGTH) seqPosition+=length;
                    else {
                        // INS is assumed to add characters to the right of
                        // $referencePosition$. We only add INS events that fall
                        // strictly inside the reference window.
						if (referencePosition>=referenceTo) break;
						else if (referencePosition<referenceFrom) {
                            seqPosition+=length;
                            continue;
                        }
                        addEvent(INS_POS,referencePosition,length,alignmentID);
                        addKmerVector(seq,seqPosition+1,seqPosition+length,alignmentID,canonizeKmers);
                        seqPosition+=length;
                    }
                }
			}
		}
        
        
		public String toString() {
			return "["+startA+".."+endA+"] x ["+startB+".."+endB+"] "+seq;
		}
	}
    
    
    private static final void addEvent(int type, int pos, int length, int alignmentID) {
        if (type==DEL_START || type==DEL_END) {
            if (del_last==del.length-1) {
                int[] newArray = new int[del.length<<1];
                System.arraycopy(del,0,newArray,0,del.length);
                del=newArray;
            }
            del[++del_last]=type;
            del[++del_last]=pos;
            del[++del_last]=length;
            del[++del_last]=alignmentID;
        }
        else if (type==INS_POS) {
            if (ins_last==ins.length-1) {
                int[] newArray = new int[ins.length<<1];
                System.arraycopy(ins,0,newArray,0,ins.length);
                ins=newArray;
            }
            ins[++ins_last]=pos;
            ins[++ins_last]=length;
            ins[++ins_last]=alignmentID;
        }
        else if (type==CLIP_LEFT || type==CLIP_RIGHT) {
            if (clip_last==clip.length-1) {
                int[] newArray = new int[clip.length<<1];
                System.arraycopy(clip,0,newArray,0,clip.length);
                clip=newArray;
            }
            clip[++clip_last]=type;
            clip[++clip_last]=pos;
            clip[++clip_last]=alignmentID;
        }
    }
    
    
	/**
	 * @param out output array: 
     * 0: length of the soft clip at the beginning of $cigar$; 
     * 1: the first position of $cigar$ from which to continue reading it after 
     *    the soft clip;
     * 2: length of the soft plus hard clip at the beginning of $cigar$.
	 */
	private static final void cigar_prefixClip(String cigar, int[] out) {
		char c;
		int i, p;
		int length, clip, softClip;
		final int cigarLength = cigar.length();
		
		p=0; i=0; clip=0; softClip=0;
		while (i<cigarLength) {
			c=cigar.charAt(i);
			if (!Character.isDigit(c)) {
                if (c=='H') {
                    clip+=Integer.parseInt(cigar.substring(p,i));
                    p=i+1;
                }
				else if (c=='S') {
                    length=Integer.parseInt(cigar.substring(p,i));
                    softClip=length; clip+=length;
                    p=i+1;
				    break;
                }
                else break;
			}
			i++;
		}
		out[0]=softClip; out[1]=p; out[2]=clip;
	}
	
	
	/**
	 * @param out output array: 
     * 0: length of the soft plus hard clip at the end of $cigar$;
     * 1: the last position of $cigar$ before the clip.
	 */
	private static final void cigar_suffixClip(String cigar, int[] out) {
		boolean found;
		char c;
		int i, iPrime;
		int clip;
		final int cigarLength = cigar.length();
	
		// First hard clip, if any.
        clip=0;
		c=cigar.charAt(cigarLength-1);
		if (c=='H') {
            // Hard clip
			i=cigarLength-2;
			while (i>=0) {
				c=cigar.charAt(i);
				if (!Character.isDigit(c)) break;
				i--;
			}
            clip+=Integer.parseInt(cigar.substring(i+1,cigarLength));
			// Following soft clip, if any.
			c=cigar.charAt(i);
			if (c=='S') {
				iPrime=i; i--;
				while (i>=0) {
					c=cigar.charAt(i);
					if (!Character.isDigit(c)) break;
					i--;
				}
				clip+=Integer.parseInt(cigar.substring(i+1,iPrime));
			}
		}
		else if (c=='S') {
			// Soft clip
			i=cigarLength-2;
			while (i>=0) {
				c=cigar.charAt(i);
				if (!Character.isDigit(c)) break;
				i--;
			}
			clip+=Integer.parseInt(cigar.substring(i+1,cigarLength-1));
		}
		else i=cigarLength-1;
		out[0]=clip; out[1]=i;
	}
	
	
	/**
	 * @param out output array: 
     * 0: the last position before any soft and hard clip at the end of $cigar$,
     *    if any (call this position $to$);
	 * 1: length of the alignment encoded in $cigar[from..to]$, projected on the
	 *    reference; 
     * 2: length of the alignment encoded in $cigar[from..to]$, projected on the
     *    read.
	 */
	private static final void cigar_length(String cigar, int from, int[] out) {
		char c;
		int i, p;
		int length, lengthA, lengthB;
		final int cigarLength = cigar.length();
		
		p=from; lengthA=0; lengthB=0;
		for (i=from; i<cigarLength; i++) {
			c=cigar.charAt(i);
			if (c=='M' || c=='X' || c=='=') {
				length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
				lengthA+=length; lengthB+=length;
			}
			else if (c=='I') {
				length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
				lengthB+=length;
			}
			else if (c=='D' || c=='N') {
				length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
				lengthA+=length;
			}
			else if (c=='S') break;
		}
		out[0]=p-1; out[1]=lengthA; out[2]=lengthB;
	}
    
    
    /**
     * @return -1 if $str$ does not represent a standard contig.
     */
    private static final int string2contig(String str) {
        char c;
        int length, out;
        
        out=-1;
        length=str.length();
        if (length>3 && str.substring(0,3).equalsIgnoreCase("chr")) { str=str.substring(3); length-=3; }
        if (length==1) {
            c=str.charAt(0);
            if (c=='X' || c=='x') return 23;
            else if (c=='Y' || c=='y') return 24;
            else if (c=='M' || c=='m') return 25;
            else return Integer.parseInt(str);
        }
        else if (length==2) {
            try { out=Integer.parseInt(str); return out; }
            catch (Exception e) { return -1; }
        }
        else return -1;
    }
    
    
    
    
    // -------------------------- KMER PROCEDURES ------------------------------
    
    /**
     * Appends the (possibly canonized) k-mer composition vector of
     * $sequence[first..last]$ to global variable $kmerVectors$.
     *
     * Remark: the procedure uses global variables $map,kmerPool,sb$.
     *
     * @param sequence assumed to be all lowercase.
     */
    private static final void addKmerVector(String sequence, int first, int last, int alignmentID, boolean canonized) {
        int i, j, k;
        int nKmers;
        String str;
        Kmer key, value;
        
        // Counting k-mers
        map.clear();
        if (last-first+1<KMER_LENGTH) return;
        j=-1;
        for (i=first; i<=last+1-KMER_LENGTH; i++) {
            j++;
            if (j==kmerPool.length) {
                Kmer[] newArray = new Kmer[kmerPool.length<<1];
                System.arraycopy(kmerPool,0,newArray,0,kmerPool.length);
                for (k=kmerPool.length; k<newArray.length; k++) newArray[k] = new Kmer();
                kmerPool=newArray;
            }
            key=kmerPool[j];
            key.sequence=canonized?canonize(sequence,i,i+KMER_LENGTH-1):sequence.substring(i,i+KMER_LENGTH);
            value=map.get(key);
            if (value==null) { key.count=1; map.put(key,key); }
            else value.count++;
        }
        
        // Building the vector
        nKmers=map.size();
        if (kmerVector.length<nKmers) kmerVector = new Kmer[nKmers];
        map.keySet().toArray(kmerVector);
        Arrays.sort(kmerVector,0,nKmers);
        addKmerVector_impl(kmerVector,nKmers,alignmentID);
    }
    
    
    private static final void addKmerVector_impl(Kmer[] vector, int nKmers, int alignmentID) {
        int i;
        
        kmerVectors_last++;
        if (kmerVectors_last==kmerVectors.length) {
            Kmer[][] newArray = new Kmer[kmerVectors.length<<1][0];
            System.arraycopy(kmerVectors,0,newArray,0,kmerVectors.length);
            kmerVectors=newArray;
        }
        kmerVectors[kmerVectors_last] = new Kmer[nKmers];
        for (i=0; i<nKmers; i++) kmerVectors[kmerVectors_last][i]=vector[i].clone();
        if (kmerVectors_last==kmerVectors_alignmentID.length) {
            int[] newArray = new int[kmerVectors_alignmentID.length<<1];
            System.arraycopy(kmerVectors_alignmentID,0,newArray,0,kmerVectors_alignmentID.length);
            kmerVectors_alignmentID=newArray;
        }
        kmerVectors_alignmentID[kmerVectors_last]=alignmentID;
    }
    
    
    /**
     * @param str assumed to be all lowercase;
     * @return the lexicographically smallest of $str[first..last]$ and its 
     * reverse-complement.
     */
    private static final String canonize(String str, int first, int last) {
        int i;
        char c, d;
        final int length = last-first+1;
        
        sb.delete(0,sb.length());
        for (i=first; i<=last; i++) sb.append(str.charAt(i));
        for (i=last; i>=first; i--) {
            c=str.charAt(i);
            switch(c) {
                case 'a': d='t'; break;
                case 'c': d='g'; break;
                case 'g': d='c'; break;
                case 't': d='a'; break;
                default: d='n'; break;
            }
            sb.append(d);
        }
        for (i=0; i<length; i++) {
            c=sb.charAt(i); d=sb.charAt(length+i);
            if (c<d) return sb.substring(0,length);
            else if (d<c) return sb.substring(length); 
        }
        return sb.substring(0,length);
    }
    
    
    /**
     * @return the Euclidean distance between the sorted composition vectors at 
     * positions $id1$ and $id2$ in $kmerVectors$.
     */
    private static final double kmerDistance(int id1, int id2) {
		int i1, i2, p;
        double out;
        final int last1 = kmerVectors[id1].length-1;
        final int last2 = kmerVectors[id2].length-1;
        final Kmer[] v1 = kmerVectors[id1];
        final Kmer[] v2 = kmerVectors[id2];
	    
        out=0.0; i1=0; i2=0;
		while (i1<=last1 && i2<=last2) {
            p=v1[i1].compareTo(v2[i2]);
			if (p<0) { out+=v1[i1].count*v1[i1].count; i1++; }
			else if (p>0) { out+=v2[i2].count*v2[i2].count; i2++; }
			else {
				out+=(v1[i1].count-v2[i2].count)*(v1[i1].count-v2[i2].count);
				i1++; i2++;
			}
		}
		while (i1<=last1) { out+=v1[i1].count*v1[i1].count; i1++; }
		while (i2<=last2) { out+=v2[i2].count*v2[i2].count; i2++; }
		return Math.sqrt(out);
    }
    
    
    private static class Kmer implements Comparable {
        public String sequence;
        public double count;
        
        public Kmer() { sequence=null; count=0; }
        
        public int hashCode() { return sequence.hashCode(); }
        
        public int compareTo(Object other) {
            Kmer otherKmer = (Kmer)other;
            return sequence.compareTo(otherKmer.sequence);
        }
        
        public boolean equals(Object other) {
            Kmer otherKmer = (Kmer)other;
            return sequence.equals(otherKmer.sequence);
        }
        
        public Kmer clone() {
            Kmer out = new Kmer();
            out.sequence=sequence; out.count=count;
            return out;
        }
        
        public String toString() { return sequence+":"+count; }
    }
    
}
