import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class MakeMultiallelic {
	
	private static final String END_STR = "END=";
	private static final int END_STR_LENGTH = END_STR.length();
	private static final String SVTYPE_STR = "SVTYPE=";
	private static final int SVTYPE_STR_LENGTH = SVTYPE_STR.length();
	private static final String SVLEN_STR = "SVLEN=";
	private static final int SVLEN_STR_LENGTH = SVLEN_STR.length();
	private static final String SEPARATOR = ";";
	private static final String INS = "INS";
	private static final String CHR_STR = "chr";
	private static final int CHR_STR_LENGTH = CHR_STR.length();
    
    /**
     * As in the Icelanders' paper.
     */
	private static final int IDENTITY_THRESHOLD = 250;
    
    /**
     * Long calls can be merged into a multiallelic record only with other long
     * calls.
     */
	private static final int MAX_SHORT_SV_LENGTH = 50000;
    
    /**
     * Reference chromosomes
     */
    private static StringBuilder[] chromosomes;


	public static void main(String[] args) throws IOException {
		final String INPUT_VCF = args[0];
		final String REFERENCE_FILE = args[1];
		final String OUTPUT_PREFIX = args[2];
		
		int i, n, p, q;
		int first, last, currentFirst_short, currentFirst_long, currentLast_short, currentLast_long, lastCall_short, lastCall_long;
		String str, type, currentChromosome_short, currentChromosome_long;
        StringBuilder sb1, sb2;
		BufferedReader br;
		BufferedWriter bw1, bw2, bw3, bw4;
		String[] tokens;
		String[][] calls_short, calls_long;
		
		System.err.println("Loading the reference...");
		chromosomes = new StringBuilder[24];
		for (i=0; i<chromosomes.length; i++) chromosomes[i] = new StringBuilder();
		br = new BufferedReader(new FileReader(REFERENCE_FILE));
		str=br.readLine(); i=-1;
		while (str!=null) {
			if (str.charAt(0)=='>') {
				i=chr2id(str)-1;
				if (i>=0) System.err.println("Loading chromosome "+str+"...");
			}
			else if (i>=0) chromosomes[i].append(str);
			str=br.readLine();
		}
		br.close();
		
		System.err.println("Making the VCF multiallelic...");
		calls_short = new String[1000][8];  // Arbitrary
		lastCall_short=-1;
		calls_long = new String[1000][8];  // Arbitrary
		lastCall_long=-1;
        sb1 = new StringBuilder(); sb2 = new StringBuilder();
		br = new BufferedReader(new FileReader(INPUT_VCF));
        bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"-header.txt"));
        bw2 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"-short.txt"));
        bw3 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"-long.txt"));
        bw4 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"-bnd.txt"));
		str=br.readLine(); n=0;
		currentChromosome_short=""; currentFirst_short=-1; currentLast_short=-1;
        currentChromosome_long=""; currentFirst_long=-1; currentLast_long=-1;
		while (str!=null) {
			if (str.charAt(0)=='#') {
                bw1.write(str); bw1.newLine();
				str=br.readLine(); n++;
				continue;
			}
			tokens=str.split("\t");
			first=Integer.parseInt(tokens[1]);
			p=tokens[7].indexOf(SVTYPE_STR);
			q=tokens[7].indexOf(SEPARATOR,p+SVTYPE_STR_LENGTH);
			if (q<0) q=tokens[7].length();
			type=tokens[7].substring(p+SVTYPE_STR_LENGTH,q);
            if (type.equalsIgnoreCase("BND")) {
                bw4.write(str); bw4.newLine();
				str=br.readLine(); n++;
				continue;
            }
            last=getEnd(tokens,type,true);
            if (last-first<=MAX_SHORT_SV_LENGTH) {
                last=getEnd(tokens,type,false);
    			if (!tokens[0].equalsIgnoreCase(currentChromosome_short) || first>currentLast_short+IDENTITY_THRESHOLD) {
                    printMultiallelic(calls_short,lastCall_short,bw2,sb1,sb2);
    				currentChromosome_short=tokens[0];
                    currentFirst_short=first; currentLast_short=last;
                    lastCall_short=-1;
    			}
    			else if (last>currentLast_short) currentLast_short=last;
    			lastCall_short++;
    			if (lastCall_short==calls_short.length) {
    				String[][] newArray = new String[calls_short.length<<1][0];
    				System.arraycopy(calls_short,0,newArray,0,calls_short.length);
    				for (i=calls_short.length; i<newArray.length; i++) newArray[i] = new String[8];
    				calls_short=newArray;
    			}
    			System.arraycopy(tokens,0,calls_short[lastCall_short],0,8);
            }
            else {
                last=getEnd(tokens,type,false);
    			if (!tokens[0].equalsIgnoreCase(currentChromosome_long) || first>currentLast_long+IDENTITY_THRESHOLD) {
                    printMultiallelic(calls_long,lastCall_long,bw3,sb1,sb2);
    				currentChromosome_long=tokens[0]; 
                    currentFirst_long=first; currentLast_long=last; 
                    lastCall_long=-1;
    			}
    			else if (last>currentLast_long) currentLast_long=last;
    			lastCall_long++;
    			if (lastCall_long==calls_long.length) {
    				String[][] newArray = new String[calls_long.length<<1][0];
    				System.arraycopy(calls_long,0,newArray,0,calls_long.length);
    				for (i=calls_long.length; i<newArray.length; i++) newArray[i] = new String[8];
    				calls_long=newArray;
    			}
    			System.arraycopy(tokens,0,calls_long[lastCall_long],0,8);
            }
            if (n%1000==0) System.err.println("Processed "+n+" VCF calls");
			str=br.readLine(); n++;
		}
		br.close();
        printMultiallelic(calls_short,lastCall_short,bw2,sb1,sb2); bw2.close();
        printMultiallelic(calls_long,lastCall_long,bw3,sb1,sb2); bw3.close();
        bw1.close(); bw4.close();
	}

    
	/**
	 * Remark: all positions in $calls$ are assumed to be one-based. 
     *
     * Remark: every DUP is converted to an INS that occurs at the first
     * position of the DUP.
     *
     * @param sb* temporary space.
	 */
	private static final void printMultiallelic(String[][] calls, int lastCall, BufferedWriter bw, StringBuilder sb1, StringBuilder sb2) throws IOException {
        if (lastCall==-1) return;
        final int chrID = chr2id(calls[0][0]);
        final int minPos = Integer.parseInt(calls[0][1]);
        final int maxPos = getMaxPos(calls,lastCall);
        final String ref = chromosomes[chrID-1].substring(minPos-1,maxPos);
        mergeAlts(calls,lastCall,ref,minPos,chrID-1,sb1,sb2);
        final String alt = sb1.toString();
        mergeInfos(calls,lastCall,sb1);
        final String info = sb1.toString();
        bw.write(calls[0][0]+"\t"+minPos+"\t"+calls[0][2]+"\t"+ref+"\t"+alt+"\t"+mergeQuals(calls,lastCall)+"\t"+mergeFilters(calls,lastCall)+"\t"+info);
        bw.newLine();
	}
    
    
    private static final int getMaxPos(String[][] calls, int lastCall) {
        int i;
        int max, end;
        
        max=0;
        for (i=0; i<=lastCall; i++) {
            end=getEnd(calls[i],null,false);
            if (end>max) max=end;
        }
        return max;
    }
    
    
    /**
     * Stores in $out$ the ALT allele that results from merging $calls[0..
     * lastCall].
     *
     * @param refPos the first position of $ref$ on its chromosome (1-based);
     * @param chrID zero-based;
     * @param sb temporary space.
     */
    private static final void mergeAlts(String[][] calls, int lastCall, String ref, int refPos, int chrID, StringBuilder out, StringBuilder sb) {
        boolean empty;
        int i, j, p, q;
        int pos, end;
        String svType;
        
        out.delete(0,out.length()); empty=true;
        for (i=0; i<=lastCall; i++) {
            p=calls[i][7].indexOf(SVTYPE_STR);
    		q=calls[i][7].indexOf(SEPARATOR,p+SVTYPE_STR_LENGTH);
            if (q<0) q=calls[i][7].length();
            svType=calls[i][7].substring(p+SVTYPE_STR_LENGTH,q);
            sb.delete(0,sb.length());
            if (svType.equalsIgnoreCase("DEL")) {
                pos=Integer.parseInt(calls[i][1]); end=getEnd(calls[i],svType,true);
                sb.append(ref.substring(0,pos+1-refPos));
                sb.append(ref.substring(end+1-refPos));
            }
            else if (svType.equalsIgnoreCase("INS")) {
                if (calls[i][4].indexOf("<")>=0) sb.append(calls[i][4]);
                else {
                    pos=Integer.parseInt(calls[i][1]);
                    sb.append(ref.substring(0,pos+1-refPos));
                    sb.append(calls[i][4]);
                    sb.append(ref.substring(pos+1-refPos));
                }
            }
            else if (svType.equalsIgnoreCase("INV")) {
                pos=Integer.parseInt(calls[i][1]); end=getEnd(calls[i],svType,true);
                sb.append(ref.substring(0,pos+1-refPos));
                for (j=end; j>=pos+1; j--) sb.append(complement(ref.charAt(j-refPos)));
                sb.append(ref.substring(end+1-refPos));
            }
            else if (svType.equalsIgnoreCase("DUP")) {
                pos=Integer.parseInt(calls[i][1]); end=getEnd(calls[i],svType,true);
                sb.append(ref.substring(0,pos+1-refPos));
                sb.append(chromosomes[chrID].substring(pos,end));
                sb.append(ref.substring(pos+1-refPos));
            }
            else sb.append(calls[i][4]);
            out.append((empty?"":",")+sb); empty=false;
        }
    }
    
    
    /**
     * @return the max over all QUAL fields.
     */
    private static final String mergeQuals(String[][] calls, int lastCall) {
        int i;
        double qual, maxQual;
        
        maxQual=-1.0;
        for (i=0; i<=lastCall; i++) {
            qual=-1.0;
            try { qual=Double.parseDouble(calls[i][5]); } catch (Exception e) { }
            if (qual>maxQual) maxQual=qual;
        }
        return maxQual==-1.0?".":(""+maxQual);
    }
    
    
    /**
     * @return PASS if at least one call passed the filters; the first filter
     * value otherwise.
     */
    private static final String mergeFilters(String[][] calls, int lastCall) {
        int i;
        String out;
        
        out=null;
        for (i=0; i<=lastCall; i++) {
            if (calls[i][6].equalsIgnoreCase("PASS") || calls[i][6].equalsIgnoreCase(".")) return "PASS";
            else if (out==null) out=calls[i][6];
        }
        return out;
    }
    
    
    /**
     * Stores in $out$ the merge of all INFO fields.
     */
    private static final void mergeInfos(String[][] calls, int lastCall, StringBuilder out) {
        boolean empty;
        int i;
        String str;
        
        out.delete(0,out.length()); empty=true;
        for (i=0; i<=lastCall; i++) {
            str=calls[i][7].replace("DUP","INS");
            str=str.replace("=",i+"=");
            out.append((empty?"":SEPARATOR)+str); empty=false;
        }
    }
    
    
    /**
     * @param svType recomputed by the procedure if NULL;
     * @param dupEndMode TRUE=the true end of the DUP; FALSE=the end of the DUP
     * when it is interpreted as an INS;
     * @return the last position of the reference used by the SV (1-based).
     */
    private static final int getEnd(String[] tokens, String svType, boolean dupEndMode) {
        int p, q;
        int length;
        
        // INS and DUP
        if (svType==null) {
            p=tokens[7].indexOf(SVTYPE_STR);
    		q=tokens[7].indexOf(SEPARATOR,p+SVTYPE_STR_LENGTH);
            if (q<0) q=tokens[7].length();
            svType=tokens[7].substring(p+SVTYPE_STR_LENGTH,q);
        }
        if (svType.equalsIgnoreCase("INS") || (svType.equalsIgnoreCase("DUP") && !dupEndMode)) return Integer.parseInt(tokens[1]);
        
        // Other types
		p=tokens[7].indexOf(END_STR);
        if (p>=0) {
            p+=END_STR_LENGTH;
    		q=tokens[7].indexOf(SEPARATOR,p);
    		if (q<0) q=tokens[7].length();
            return Integer.parseInt(tokens[7].substring(p,q));
        }
        else {
            p=tokens[7].indexOf(SVLEN_STR)+SVLEN_STR_LENGTH;
    		q=tokens[7].indexOf(SEPARATOR,p);
    		if (q<0) q=tokens[7].length();
            length=Integer.parseInt(tokens[7].substring(p,q));
            if (length<0) length=-length;
            return Integer.parseInt(tokens[1])+length;
        }
    }
    
    
	/**
	 * @return one-based (zero=all normal chromosomes have been read).
	 */
	private static final int chr2id(String str) {
		final int p = str.indexOf(CHR_STR);
		int q = str.indexOf(" ",p+CHR_STR_LENGTH+1);
		if (q<0) q=str.length();
		if (q==p+CHR_STR_LENGTH+1) {
			final char c = str.charAt(p+CHR_STR_LENGTH);
			if (c=='X' || c=='x') return 23;
			else if (c=='Y' || c=='y') return 24;
			else if (c=='M' || c=='m') return 0;
			else return Integer.parseInt(c+"");
		}
		else {	
			int out = 0;
			try { out=Integer.parseInt(str.substring(p+CHR_STR_LENGTH,q)); }
			catch (NumberFormatException e) { }
			return out;
		}
	}
    
    
    private static final char complement(char c) {
        switch (c) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'a': return 't';
            case 'c': return 'g';
            case 'g': return 'c';
            case 't': return 'a';
            default: return 'N';
        }
    }
	
}