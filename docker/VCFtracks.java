import java.util.Arrays;
import java.util.HashMap;
import java.io.*;


/**
 * 
 */
public class VCFtracks {
	/**
	 * SVs shorter than this are ignored
	 */
	private static final int MIN_SV_LENGTH = 40;
    private static final int INS_SLACK = MIN_SV_LENGTH>>1;  // Arbitrary
	private static final boolean ONLY_PASS = true;
	
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
     * SV types that produce events
	 */
	public static final byte TYPE_DELETION = 0;
    public static final byte TYPE_INVERSION = 1;
    public static final byte TYPE_DUPLICATION = 2;
    public static final byte TYPE_INSERTION = 3;
    
    /**
     * Event types
     */
    private static final int DEL_START = 0;
    private static final int DEL_END = 1;
    private static final int INV_START = 2;
    private static final int INV_END = 3;
    private static final int DUP_START = 4;
    private static final int DUP_END = 5;
    private static final int INS_POS = 6;
    
    /**
     * Number of signals reported in output
     */
    public static final int N_SIGNALS = 16;
	
	/**
	 * All the VCF records that intersect the current window
	 */
	private static Call[] calls;
	private static int lastCall;
    
    /**
     * All events in the current window
     */
    private static int[] del, inv, dup, ins;
    private static int del_last, inv_last, dup_last, ins_last;
    private static Kmer[][] kmerVectors;
    private static int[] kmerVectors_callID;
    private static int kmerVectors_last;
	
	/**
	 * Temporary space
	 */
    private static StringBuilder sb;
    private static HashMap<Kmer,Kmer> map;
    private static Kmer[] kmerPool, kmerVector;
    private static int[] coverageHistogram;
	
	
	/**
	 * Given a region of a chromosome, prints window-based statistics over it.
     * Should be used on the merged VCF of a cohort.
     *
     * GRCh38:
     * chr21 length=46709983
     * chr22 length=50818468
     *
     * T2T:
     * chr21 length=45090682
     * chr22 length=51324926
	 */
	public static void main(String[] args) throws IOException {
		int i, p, q;
		int contig, position, type, length, currentContig, currentStart, nCalls;
		String str, alt, filter, info, tmpString;
		BufferedReader br;
		BufferedWriter bw;
		
		// Parsing the input
		final String VCF_FILE = args[0];  // Sorted
        final String CHROMOSOME_ID = args[1];
        final int FROM_POS = Integer.parseInt(args[2]);  // 1-based
        final int TO_POS = Integer.parseInt(args[3]);  // 1-based
        WINDOW_LENGTH=Integer.parseInt(args[4]);
        WINDOW_STEP=Integer.parseInt(args[5]);
        KMER_LENGTH=Integer.parseInt(args[6]);
		final String OUTPUT_FILE = args[7];
		
		// Allocating memory
        currentContig=string2contig(CHROMOSOME_ID);
		calls = new Call[50];  // Arbitrary
		for (i=0; i<calls.length; i++) calls[i] = new Call();
        del = new int[300];  // Arbitrary
        inv = new int[300];  // Arbitrary
        dup = new int[300];  // Arbitrary
        ins = new int[200];  // Arbitrary
        kmerVectors = new Kmer[100][0];  // Arbitrary
        kmerVectors_callID = new int[100];  // Arbitrary
        sb = new StringBuilder();
        map = new HashMap<Kmer,Kmer>();
        kmerPool = new Kmer[100];  // Arbitrary
        for (i=0; i<kmerPool.length; i++) kmerPool[i] = new Kmer();
        kmerVector = new Kmer[100];  // Arbitrary
        coverageHistogram = new int[WINDOW_LENGTH];
		
		// Scanning the reference
		br = new BufferedReader(new FileReader(VCF_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		nCalls=0; currentStart=FROM_POS; lastCall=-1; 
		str=br.readLine();
		while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
			nCalls++;
			p=str.indexOf(VCF_SEPARATOR); 
            contig=string2contig(str.substring(0,p));
            q=str.indexOf(VCF_SEPARATOR,p+1);
            position=Integer.parseInt(str.substring(p+1,q));  // 1-based
			if (contig<currentContig) {
                str=br.readLine();
                continue;
			}
            else if (contig>currentContig || position>TO_POS) {
                getTracks(currentContig,currentStart,bw,CANONIZE_KMERS);
                currentStart+=WINDOW_STEP;
                break;
            }
            while (position>currentStart+WINDOW_LENGTH-1) {
				getTracks(currentContig,currentStart,bw,CANONIZE_KMERS);
				currentStart+=WINDOW_STEP;
				filterCalls(currentStart);
            }
            p=str.indexOf(VCF_SEPARATOR,q+1);
            p=str.indexOf(VCF_SEPARATOR,p+1); q=str.indexOf(VCF_SEPARATOR,p+1);
            alt=str.substring(p+1,q);
            p=str.indexOf(VCF_SEPARATOR,q+1); q=str.indexOf(VCF_SEPARATOR,p+1);
            filter=str.substring(p+1,q);
            p=q; q=str.indexOf(VCF_SEPARATOR,p+1);
            info=str.substring(p+1,q);
            tmpString=getField(info,SVTYPE_STR);
            
// System.err.println("type="+tmpString+" info field="+info);
// if (tmpString!=null) System.err.println("which is number "+getType_infoField(tmpString));
            
            if (tmpString!=null) type=getType_infoField(tmpString);
            else if (alt.charAt(0)=='<') type=getType_infoField(alt.substring(1,alt.length()-1));
            else type=-1;
            tmpString=getField(info,SVLEN_STR);
            if (tmpString!=null) {
                length=Integer.parseInt(tmpString);
                if (length<0) length=-length;
            }
            else if (type==TYPE_INSERTION && alt.charAt(0)!='<') length=alt.length()-1;
            else {
                tmpString=getField(info,END_STR);
                if (tmpString!=null) length=Integer.parseInt(tmpString)-(position+1)+1;
                else length=-1;
            }
			appendCall(position,alt,filter,type,length,FROM_POS,TO_POS);
			if (nCalls%10000==0) System.err.println("Processed "+nCalls+" calls");
			str=br.readLine();
		}
        while (currentStart<=TO_POS) {
            bw.write(currentContig+","+currentStart);
            for (i=0; i<N_SIGNALS; i++) bw.write(",0");
            bw.newLine();
            currentStart+=WINDOW_STEP;
        }
		br.close(); bw.close();
	}
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	public static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}
    
    
	/**
	 * @return -1 iff the type cannot be determined.
	 */
	private static final int getType_infoField(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return TYPE_DELETION;
		else if ( type.equalsIgnoreCase(INS_STR) || 
			      type.equalsIgnoreCase(INS_ME_STR) || 
				  type.equalsIgnoreCase(INS_NOVEL_STR)
				) return TYPE_INSERTION;
		else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return TYPE_DUPLICATION;
		else if (type.equalsIgnoreCase(INV_STR)) return TYPE_INVERSION;
		else return -1;
	}
    
    
    /**
     * VCF constants
     */
    public static final char COMMENT = '#';
    private static final String VCF_SEPARATOR = "\t";
    private static final String PASS_STR = "PASS";
	public static final String END_STR = "END";
    public static final String CIEND_STR = "CIEND";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
    public static final String SEPARATOR = ";";
    
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
	public static final String INS_NOVEL_STR = "INS:NOVEL";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";

	
    /**
     * @param fromPos,toPos calls that do not intersect $[fromPos..toPos]$ in 
     * the reference (1-based) are not appended;
     * @param type if not supported, the call is discarded;
     * @param length if < MIN_SV_LENGTH, the call is discarded.
     */
	private static final void appendCall(int position, String alt, String filter, int type, int length, int fromPos, int toPos) {
		int i;
        
// if (type==TYPE_INVERSION || type==TYPE_DUPLICATION) System.err.println("appendCall> position="+position+" alt="+alt+" filter="+filter+" type="+type+" length="+length+" fromPos="+fromPos+" toPos="+toPos);
        
		if ((ONLY_PASS && !filter.equalsIgnoreCase(PASS_STR) && !filter.equals(".")) || type==-1 || length<MIN_SV_LENGTH) return;
		lastCall++;
		if (lastCall==calls.length) {
			Call[] newArray = new Call[calls.length<<1];
			System.arraycopy(calls,0,newArray,0,calls.length);
			for (i=calls.length; i<newArray.length; i++) newArray[i] = new Call();
			calls=newArray;
		}
		calls[lastCall].set(position,alt,filter,type,length);
        if (calls[lastCall].end<fromPos || calls[lastCall].start>toPos) lastCall--;
	}
	
	
	/**
	 * Removes from $calls$ every call whose projection on the reference does 
     * not intersect $[from..from+WINDOW_LENGTH-1]$ (1-based).
	 */
	private static final void filterCalls(int from) {
		int i, j;
		final int to = from+WINDOW_LENGTH-1;
		Call tmpCall;
		
		j=-1;
		for (i=0; i<=lastCall; i++) {
			if (Math.min(to,calls[i].end)>=Math.max(from,calls[i].start)) {
				j++;
				tmpCall=calls[i];
				calls[i]=calls[j];
				calls[j]=tmpCall;
			}
		}
		lastCall=j;
	}
    
    
    /**
     * Given window $currentContig[currentStart..currentStart+WINDOW_LENGTH-1]$,
     * for each event type, the procedure outputs the max distance between an 
     * event in the window, and its nearest neighbor in the window. The 
     * procedure outputs 0 iff there is <=1 event of a given type in the window.
     *
     * @param currentStart 1-based.
     */
    private static final void getTracks(int currentContig, int currentStart, BufferedWriter bw, boolean canonizeKmers) throws IOException {
        int i, j;
        int dx, dy;
        int nDel, nDelStart, nDelEnd, nInv, nInvStart, nInvEnd, nDup, nDupStart, nDupEnd, nIns;
        double sum, distance, minDistance, maxDel, maxInv, maxDup, maxIns, maxKmer, avgCoverage;
        

// System.err.println("getTracks> Calls in ["+currentStart+"..:");
// for (int x=0; x<=lastCall; x++) System.err.println(calls[x]);
        
        
        if (lastCall==-1) {
            bw.write(currentContig+","+currentStart);
            for (i=0; i<N_SIGNALS; i++) bw.write(",0");
            bw.newLine();
            return;
        }
                
        // Collecting events, in no particular order.
        nDel=0; nInv=0; nDup=0;
        del_last=-1; inv_last=-1; dup_last=-1; ins_last=-1; kmerVectors_last=-1;
        Arrays.fill(coverageHistogram,0);
        for (i=0; i<=lastCall; i++) {
            calls[i].addEvents(currentStart,i,canonizeKmers,coverageHistogram);
            switch (calls[i].type) {
                case TYPE_DELETION: nDel++; break;
                case TYPE_INVERSION: nInv++; break;
                case TYPE_DUPLICATION: nDup++; break;
            }
        }
        avgCoverage=0.0;
        for (i=0; i<WINDOW_LENGTH; i++) avgCoverage+=coverageHistogram[i];
        avgCoverage/=WINDOW_LENGTH;
        
        // DEL
        maxDel=0.0; nDelStart=0; nDelEnd=0;
        for (i=0; i<=del_last; i+=3) {
            if (del[i]==DEL_START) nDelStart++;
            else if (del[i]==DEL_END) nDelEnd++;
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=del_last; j+=3) {
                if (j==i || del[j]!=del[i]) continue;
                dx=del[j+1]-del[i+1]; dy=del[j+2]-del[i+2];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  
                // Isolated event: NOP.
            }
            else if (minDistance>maxDel) maxDel=minDistance;
        }
        // INV
        maxInv=0.0; nInvStart=0; nInvEnd=0;
        for (i=0; i<=inv_last; i+=3) {
            if (inv[i]==INV_START) nInvStart++;
            else if (inv[i]==INV_END) nInvEnd++;
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=inv_last; j+=3) {
                if (j==i || inv[j]!=inv[i]) continue;
                dx=inv[j+1]-inv[i+1]; dy=inv[j+2]-inv[i+2];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  
                // Isolated event: NOP.
            }
            else if (minDistance>maxInv) maxInv=minDistance;
        }
        // DUP
        maxDup=0.0; nDupStart=0; nDupEnd=0;
        for (i=0; i<=dup_last; i+=3) {
            if (dup[i]==DUP_START) nDupStart++;
            else if (dup[i]==DUP_END) nDupEnd++;
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=dup_last; j+=3) {
                if (j==i || dup[j]!=dup[i]) continue;
                dx=dup[j+1]-dup[i+1]; dy=dup[j+2]-dup[i+2];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  
                // Isolated event: NOP.
            }
            else if (minDistance>maxDup) maxDup=minDistance;
        }
        // INS
        maxIns=0.0; nIns=ins_last+1;
        for (i=0; i<=ins_last; i+=2) {
            minDistance=Integer.MAX_VALUE;
            for (j=0; j<=ins_last; j+=2) {
                if (j==i) continue;
                dx=ins[j]-ins[i]; dy=ins[j+1]-ins[i+1];
                distance=Math.sqrt(dx*dx+dy*dy);
                if (distance<minDistance) minDistance=distance;
            }
            if (minDistance==Integer.MAX_VALUE) {  
                // Isolated event: NOP.
            }
            else if (minDistance>maxIns) maxIns=minDistance;
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
                    if (j==i) continue;
                    distance=kmerDistance(i,j); 
                    if (distance<minDistance) minDistance=distance; 
                } 
                if (minDistance==Integer.MAX_VALUE) {
                    // Isolated event: NOP.
                }
                else if (minDistance>maxKmer) maxKmer=minDistance;
            }
        }
        
        bw.write(currentContig+","+currentStart+","+maxDel+","+maxInv+","+maxDup+","+maxIns+","+maxKmer+","+avgCoverage+","+nDel+","+nDelStart+","+nDelEnd+","+nInv+","+nInvStart+","+nInvEnd+","+nDup+","+nDupStart+","+nDupEnd+","+nIns+"\n");
    }
    
    
    

	private static class Call {
		public int start, end;  // 1-based
        public String alt;  // Might be symbolic
        public String filter;
		public int type, length;
		
		
		public final void set(int position, String alt, String filter, int type, int length) {
            if (type==TYPE_INSERTION) { 
                // An INS is assigned to the position immediately to its left.
                start=position; end=position;
                if (alt.charAt(0)!='<') this.length=alt.length()-1;
                else this.length=length;
            }
			else { 
                start=position+1; end=position+length; this.length=length;
            }
            this.alt=alt; this.filter=filter; this.type=type;
		}
		
        
		/**
         * Adds to the global lists all events that fall in the window
         * $[referenceFrom .. referenceFrom+WINDOW_LENGTH-1]$.
         *
         * @param referenceFrom 1-based.
		 */
		public final void addEvents(int referenceFrom, int alignmentID, boolean canonizeKmers, int[] coverageHistogram) {
            boolean addCoverage;
            int i;
            int from, to;
			final int referenceTo = referenceFrom+WINDOW_LENGTH-1;
            
            addCoverage=false; from=-1; to=-1;
            if (type==TYPE_DELETION) {
                if (start>=referenceFrom && start<=referenceTo) addEvent(DEL_START,start,length);
                if (end>=referenceFrom && end<=referenceTo) addEvent(DEL_END,end,length);
                from=start; to=end; addCoverage=true;
            }
            else if (type==TYPE_INVERSION) {
                if (start>=referenceFrom && start<=referenceTo) addEvent(INV_START,start,length);
                if (end>=referenceFrom && end<=referenceTo) addEvent(INV_END,end,length);
                from=start; to=end; addCoverage=true;
            }
            else if (type==TYPE_DUPLICATION) {
                if (start>=referenceFrom && start<=referenceTo) addEvent(DUP_START,start,length);
                if (end>=referenceFrom && end<=referenceTo) addEvent(DUP_END,end,length);
                from=start; to=end; addCoverage=true;
            }
            else if (type==TYPE_INSERTION) {
                if (start>=referenceFrom && start<=referenceTo) addEvent(INS_POS,start,length);
                if (alt.charAt(0)!='<') addKmerVector(alt,1,length-1,canonizeKmers);
                from=start-INS_SLACK; to=start+INS_SLACK; addCoverage=true;
            }
            if (addCoverage) {
                if (from<referenceFrom) from=referenceFrom;
                if (to>referenceTo) to=referenceTo;
                for (i=from; i<=to; i++) coverageHistogram[i-referenceFrom]++;
            }
		}
        
        
		public String toString() {
			return "T="+type+" ["+start+".."+end+"] L="+length+" ALT="+alt;
		}
	}
    
    
    private static final void addEvent(int type, int pos, int length) {
        if (type==DEL_START || type==DEL_END) {
            if (del_last==del.length-1) {
                int[] newArray = new int[del.length<<1];
                System.arraycopy(del,0,newArray,0,del.length);
                del=newArray;
            }
            del[++del_last]=type;
            del[++del_last]=pos;
            del[++del_last]=length;
        }
        else if (type==INV_START || type==INV_END) {
            if (inv_last==inv.length-1) {
                int[] newArray = new int[inv.length<<1];
                System.arraycopy(inv,0,newArray,0,inv.length);
                inv=newArray;
            }
            inv[++inv_last]=type;
            inv[++inv_last]=pos;
            inv[++inv_last]=length;
        }
        else if (type==DUP_START || type==DUP_END) {
            if (dup_last==dup.length-1) {
                int[] newArray = new int[dup.length<<1];
                System.arraycopy(dup,0,newArray,0,dup.length);
                dup=newArray;
            }
            dup[++dup_last]=type;
            dup[++dup_last]=pos;
            dup[++dup_last]=length;
        }
        else if (type==INS_POS) {
            if (ins_last==ins.length-1) {
                int[] newArray = new int[ins.length<<1];
                System.arraycopy(ins,0,newArray,0,ins.length);
                ins=newArray;
            }
            ins[++ins_last]=pos;
            ins[++ins_last]=length;
        }
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
    private static final void addKmerVector(String sequence, int first, int last, boolean canonized) {
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
        addKmerVector_impl(kmerVector,nKmers);
    }
    
    
    private static final void addKmerVector_impl(Kmer[] vector, int nKmers) {
        int i;
        
        kmerVectors_last++;
        if (kmerVectors_last==kmerVectors.length) {
            Kmer[][] newArray = new Kmer[kmerVectors.length<<1][0];
            System.arraycopy(kmerVectors,0,newArray,0,kmerVectors.length);
            kmerVectors=newArray;
        }
        kmerVectors[kmerVectors_last] = new Kmer[nKmers];
        for (i=0; i<nKmers; i++) kmerVectors[kmerVectors_last][i]=vector[i].clone();
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
