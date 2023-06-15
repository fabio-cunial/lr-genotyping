import java.io.*;


/** 
 * Prints histograms of the number of windows with a given value of anomalous 
 * VCF signals (in each window, calls are merged across all samples).
 */
public class VCFtracksHistograms {
    /**
     * Histogram bins
     */
    private static final double MAXDEL_QUANTUM = 10.0;
    private static final double MAXDEL_MAX = 1000.0;
    private static final double MAXINS_QUANTUM = 10.0;
    private static final double MAXINS_MAX = 5000.0;
    private static final double MAXKMER_QUANTUM = 0.01;
    private static final double MAXKMER_MAX = 10.0;
    private static final double COV_QUANTUM = 0.5;
    private static final double COV_MAX = 200.0;
    private static final double ENTROPY_QUANTUM = 0.05;
    private static final double ENTROPY_MAX = 10.0;
    
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final String OUTPUT_PREFIX = args[1];
        
        final int N_COORDINATE_COLUMNS = 2;
        
        int i;
        int length, nAnomalous;
        String str;
        BufferedReader br;
        long[] maxDel, maxIns, maxKmer, coverageDel, coverageIns, coverageAll, entropyDel, entropyIns, entropyAll;
        String[] tokens;
        
        maxDel = new long[1+(int)Math.ceil(MAXDEL_MAX/MAXDEL_QUANTUM)];
        maxIns = new long[1+(int)Math.ceil(MAXINS_MAX/MAXINS_QUANTUM)];
        maxKmer = new long[1+(int)Math.ceil(MAXKMER_MAX/MAXKMER_QUANTUM)];
        coverageDel = new long[1+(int)Math.ceil(COV_MAX/COV_QUANTUM)];
        coverageIns = new long[1+(int)Math.ceil(COV_MAX/COV_QUANTUM)];
        coverageAll = new long[1+(int)Math.ceil(COV_MAX/COV_QUANTUM)];
        entropyDel = new long[1+(int)Math.ceil(ENTROPY_MAX/ENTROPY_QUANTUM)];
        entropyIns = new long[1+(int)Math.ceil(ENTROPY_MAX/ENTROPY_QUANTUM)];
        entropyAll = new long[1+(int)Math.ceil(ENTROPY_MAX/ENTROPY_QUANTUM)];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            incrementHistogram(Double.parseDouble(tokens[2]),maxDel,MAXDEL_MAX,MAXDEL_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[5]),maxIns,MAXINS_MAX,MAXINS_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[6]),maxKmer,MAXKMER_MAX,MAXKMER_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[7]),coverageDel,COV_MAX,COV_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[8]),coverageIns,COV_MAX,COV_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[9]),coverageAll,COV_MAX,COV_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[10]),entropyDel,ENTROPY_MAX,ENTROPY_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[11]),entropyIns,ENTROPY_MAX,ENTROPY_QUANTUM);
            incrementHistogram(Double.parseDouble(tokens[12]),entropyAll,ENTROPY_MAX,ENTROPY_QUANTUM);
            str=br.readLine();
        }
        br.close();
        printHistogram(maxDel,MAXDEL_QUANTUM,OUTPUT_PREFIX+"-maxDel.txt");
        printHistogram(maxIns,MAXINS_QUANTUM,OUTPUT_PREFIX+"-maxIns.txt");
        printHistogram(maxKmer,MAXKMER_QUANTUM,OUTPUT_PREFIX+"-maxKmer.txt");
        printHistogram(coverageDel,COV_QUANTUM,OUTPUT_PREFIX+"-coverage-del.txt");
        printHistogram(coverageIns,COV_QUANTUM,OUTPUT_PREFIX+"-coverage-ins.txt");
        printHistogram(coverageAll,COV_QUANTUM,OUTPUT_PREFIX+"-coverage-all.txt");
        printHistogram(entropyDel,ENTROPY_QUANTUM,OUTPUT_PREFIX+"-entropy-del.txt");
        printHistogram(entropyIns,ENTROPY_QUANTUM,OUTPUT_PREFIX+"-entropy-ins.txt");
        printHistogram(entropyAll,ENTROPY_QUANTUM,OUTPUT_PREFIX+"-entropy-all.txt");
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

}