import java.io.*;


/** 
 * 
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
     * @param args 1 take windows only from the first $X$ samples.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int N_SAMPLES_OUTPUT = Integer.parseInt(args[1]);
        final int N_SAMPLES_TOTAL = Integer.parseInt(args[2]);
        final String OUTPUT_DIR = args[3];
        
        final int N_COORDINATE_COLUMNS = 2;
        final int N_COLUMNS = N_COORDINATE_COLUMNS+N_SAMPLES_TOTAL*BAMtracks.N_SIGNALS;
        
        int i;
        int length, nSamples;
        String str;
        BufferedReader br;
        int[] maxDel, maxIns, maxClip, maxKmer, coverage, quality;
        String[] tokens;
        
        maxDel = new int[1+(int)Math.ceil(MAXDEL_MAX/MAXDEL_QUANTUM)];
        maxIns = new int[1+(int)Math.ceil(MAXINS_MAX/MAXINS_QUANTUM)];
        maxClip = new int[1+(int)Math.ceil(MAXCLIP_MAX/MAXCLIP_QUANTUM)];
        maxKmer = new int[1+(int)Math.ceil(MAXKMER_MAX/MAXKMER_QUANTUM)];
        coverage = new int[1+(int)Math.ceil(COV_MAX/COV_QUANTUM)];
        quality = new int[1+(int)Math.ceil(QUAL_MAX/QUAL_QUANTUM)];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            length=tokens.length;
            if (length!=N_COLUMNS) {
                System.err.println("nTokens="+length);
                str=br.readLine();
                continue;
            }
            nSamples=0;
            for (i=N_COORDINATE_COLUMNS; i<N_COLUMNS; i+=BAMtracks.N_SIGNALS) {
                incrementHistogram(Double.parseDouble(tokens[i]),maxDel,MAXDEL_MAX,MAXDEL_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+1]),maxIns,MAXINS_MAX,MAXINS_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+2]),maxClip,MAXCLIP_MAX,MAXCLIP_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+3]),maxKmer,MAXKMER_MAX,MAXKMER_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+4]),coverage,COV_MAX,COV_QUANTUM);
                incrementHistogram(Double.parseDouble(tokens[i+5]),quality,QUAL_MAX,QUAL_QUANTUM);
                nSamples++;
                if (nSamples==N_SAMPLES_OUTPUT) break;
            }
            str=br.readLine();
        }
        br.close();
        printHistogram(maxDel,MAXDEL_QUANTUM,OUTPUT_DIR+"/histogram-maxDel.txt");
        printHistogram(maxIns,MAXINS_QUANTUM,OUTPUT_DIR+"/histogram-maxIns.txt");
        printHistogram(maxClip,MAXCLIP_QUANTUM,OUTPUT_DIR+"/histogram-maxClip.txt");
        printHistogram(maxKmer,MAXKMER_QUANTUM,OUTPUT_DIR+"/histogram-maxKmer.txt");
        printHistogram(coverage,COV_QUANTUM,OUTPUT_DIR+"/histogram-coverage.txt");
        printHistogram(quality,QUAL_QUANTUM,OUTPUT_DIR+"/histogram-quality.txt");
    }
    
    
    private static final void incrementHistogram(double value, int[] histogram, double max, double quantum) {
        if (value==0.0) histogram[0]++;
        else if (value>=max) histogram[histogram.length-1]++;
        else if (value>0) histogram[1+(int)(value/quantum)]++;
    }
    
    
    private static final void printHistogram(int[] histogram, double quantum, String path) throws IOException {
        int i;
        final int length = histogram.length;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        bw.write("0.0,"+histogram[0]+"\n");
        for (i=1; i<length; i++) bw.write((quantum*i)+","+histogram[i]+"\n");
        bw.close();
    }

}