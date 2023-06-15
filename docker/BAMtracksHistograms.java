import java.io.*;


/** 
 * Prints histograms of the number of windows with a given value of anomalous 
 * BAM signals, where a window is an interval of the reference in a specific
 * sample (i.e. the same interval of the reference gives rise to nSamples
 * windows). Prints also a BED file with every window that is anomalous in
 * some sample.
 */
public class BAMtracksHistograms {
    /**
     * Histogram bins
     */
    private static final double MAXDEL_QUANTUM = 10.0;
    private static final double MAXDEL_MAX = 1000.0;
    private static final double MAXINS_QUANTUM = 10.0;
    private static final double MAXINS_MAX = 5000.0;        
    private static final double MAXCLIP_QUANTUM = 10.0;
    private static final double MAXCLIP_MAX = 1000.0;
    private static final double MAXKMER_QUANTUM = 0.01;
    private static final double MAXKMER_MAX = 10.0;
    private static final double COV_QUANTUM = 0.5;
    private static final double COV_MAX = 200.0;
    private static final double QUAL_QUANTUM = 1.0;
    private static final double QUAL_MAX = 60.0;
    
    
    /**
     * @param args 
     * 1: take windows only from the first $X$ samples;
     * 9: assumed to have a one-line header.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int N_SAMPLES_OUTPUT = Integer.parseInt(args[1]);
        final int N_SAMPLES_TOTAL = Integer.parseInt(args[2]);
        final String OUTPUT_PREFIX = args[3];
        final double ANOMALOUS_MAXDEL = Double.parseDouble(args[4]);
        final double ANOMALOUS_MAXINS = Double.parseDouble(args[5]);
        final double ANOMALOUS_MAXCLIP = Double.parseDouble(args[6]);
        final double ANOMALOUS_MAXKMER = Double.parseDouble(args[7]);
        final double ANOMALOUS_COV = Double.parseDouble(args[8]);
        final int WINDOW_LENGTH = Integer.parseInt(args[9]);
        
        final int N_COORDINATE_COLUMNS = 2;
        final int N_COLUMNS = N_COORDINATE_COLUMNS+N_SAMPLES_TOTAL*BAMtracks.N_SIGNALS;
        
        int i, p;
        int length, nSamples, nAnomalous;
        String str, strPrime;
        BufferedReader br;
        BufferedWriter bw;
        int[] bed_chr;
        long[] maxDel, maxIns, maxClip, maxKmer, coverage, quality, anomalous;
        int[][] bed_intervals;
        String[] tokens;
        
        maxDel = new long[1+(int)Math.ceil(MAXDEL_MAX/MAXDEL_QUANTUM)];
        maxIns = new long[1+(int)Math.ceil(MAXINS_MAX/MAXINS_QUANTUM)];
        maxClip = new long[1+(int)Math.ceil(MAXCLIP_MAX/MAXCLIP_QUANTUM)];
        maxKmer = new long[1+(int)Math.ceil(MAXKMER_MAX/MAXKMER_QUANTUM)];
        coverage = new long[1+(int)Math.ceil(COV_MAX/COV_QUANTUM)];
        quality = new long[1+(int)Math.ceil(QUAL_MAX/QUAL_QUANTUM)];
        anomalous = new long[1+N_SAMPLES_TOTAL];
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"-anomalous.bed"));
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine(); p=0;
        while (str!=null) {
            tokens=str.split(",");
            length=tokens.length;
            if (length!=N_COLUMNS) {
                System.err.println("nTokens="+length);
                str=br.readLine();
                continue;
            }
            nSamples=0; nAnomalous=0;
            for (i=N_COORDINATE_COLUMNS; i<N_COLUMNS; i+=BAMtracks.N_SIGNALS) {
                incrementHistogram(Double.parseDouble(tokens[i]),maxDel,MAXDEL_MAX,MAXDEL_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+1]),maxIns,MAXINS_MAX,MAXINS_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+2]),maxClip,MAXCLIP_MAX,MAXCLIP_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+3]),maxKmer,MAXKMER_MAX,MAXKMER_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+4]),coverage,COV_MAX,COV_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+5]),quality,QUAL_MAX,QUAL_QUANTUM);
                if ( Double.parseDouble(tokens[i])>=ANOMALOUS_MAXDEL ||
                     Double.parseDouble(tokens[i+1])>=ANOMALOUS_MAXINS ||
                     Double.parseDouble(tokens[i+2])>=ANOMALOUS_MAXCLIP ||
                     Double.parseDouble(tokens[i+3])>=ANOMALOUS_MAXKMER ||
                     Double.parseDouble(tokens[i+4])>=ANOMALOUS_COV
                   ) nAnomalous++;
                nSamples++;
                if (nSamples==N_SAMPLES_OUTPUT) break;
            }
            anomalous[nAnomalous]++;
            if (nAnomalous>=0.2*nSamples) bw.write("chr"+tokens[0]+"\t"+tokens[1]+"\t"+(Integer.parseInt(tokens[1])+WINDOW_LENGTH-1)+"\n");
            str=br.readLine();
        }
        br.close(); bw.close();
        printHistogram(maxDel,MAXDEL_QUANTUM,OUTPUT_PREFIX+"-maxDel.txt");
        printHistogram(maxIns,MAXINS_QUANTUM,OUTPUT_PREFIX+"-maxIns.txt");
        printHistogram(maxClip,MAXCLIP_QUANTUM,OUTPUT_PREFIX+"-maxClip.txt");
        printHistogram(maxKmer,MAXKMER_QUANTUM,OUTPUT_PREFIX+"-maxKmer.txt");
        printHistogram(coverage,COV_QUANTUM,OUTPUT_PREFIX+"-coverage.txt");
        printHistogram(quality,QUAL_QUANTUM,OUTPUT_PREFIX+"-quality.txt");
        printHistogramPrime(anomalous,OUTPUT_PREFIX+"-anomalous.txt");
    }
    
    
    private static final void incrementHistogram(double value, long[] histogram, double max, double quantum) {
        if (value==0.0) histogram[0]++;
        else if (value>=max) histogram[histogram.length-1]++;
        else if (value>0) histogram[1+(int)(value/quantum)]++;
    }
    
    
    private static final void printHistogram(long[] histogram, double quantum, String path) throws IOException {
        int i;
        final int length = histogram.length;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        bw.write("0.0,"+histogram[0]+"\n");
        for (i=1; i<length; i++) bw.write((quantum*i)+","+histogram[i]+"\n");
        bw.close();
    }
    
    
    private static final void printHistogramPrime(long[] histogram, String path) throws IOException {
        int i;
        final int length = histogram.length;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        for (i=0; i<length; i++) bw.write(i+","+histogram[i]+"\n");
        bw.close();
    }
    
    
	/**
	 * @return the integer ID (one-based) of the contig that $str$ encodes, or 
	 * -1 if the contig could not be found.
	 *
	 * Contig IDs have the following order:
	 * 1 .. 22 : chromosomes;
	 * 23 : X chromosome;
	 * 24 : Y chromosome;
	 * 25 : M chromosome.
	 */
	public static final int string2contig(String str) {
		final int STR_LENGTH = str.length();
		boolean found = false;
		int i;
		
		for (i=0; i<STR_LENGTH; i++) {
			if (!Character.isDigit(str.charAt(0))) {
				found=true;
				break;
			}
		}
		if (!found) return Integer.parseInt(str);
		else if (str.equalsIgnoreCase(X_STR) || str.equalsIgnoreCase(X_STR_PRIME)) return 23;
		else if (str.equalsIgnoreCase(Y_STR) || str.equalsIgnoreCase(Y_STR_PRIME)) return 24;
		else if (str.equalsIgnoreCase(M_STR) || str.equalsIgnoreCase(M_STR_PRIME) || str.equalsIgnoreCase(MT_STR) || str.equalsIgnoreCase(MT_STR_PRIME)) return 25;
		else if (str.substring(0,CHR_STR_LENGTH).equalsIgnoreCase(CHR_STR) && str.length()<=CHR_STR_LENGTH+2) return Integer.parseInt(str.substring(CHR_STR_LENGTH));
		else return -1;
	}
    
    
	/**
	 * Chromosomes
	 */
	public static final int N_CHROMOSOMES = 25;  // 22+X+Y+M
	public static final String CHR_STR = "chr";
	public static final String X_STR_PRIME = "X";
	public static final String Y_STR_PRIME = "Y";
	public static final String M_STR_PRIME = "M";
	public static final String MT_STR_PRIME = "MT";
	public static final int CHR_STR_LENGTH = CHR_STR.length();
	public static final String X_STR = CHR_STR+X_STR_PRIME;
	public static final String Y_STR = CHR_STR+Y_STR_PRIME;
	public static final String M_STR = CHR_STR+M_STR_PRIME;
	public static final String MT_STR = CHR_STR+MT_STR_PRIME;

}