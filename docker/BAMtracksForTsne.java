import java.io.*;
import java.awt.*;
import java.awt.image.*;
import javax.imageio.*;


/** 
 * Prints in output the union of all windows with an event (not subdivided by 
 * sample).
 */
public class BAMtracksForTsne {
    
    /**
     * @param args 2 take windows only from the first $X$ samples.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int N_SAMPLES_OUTPUT = Integer.parseInt(args[1]);
        final int N_SAMPLES_TOTAL = Integer.parseInt(args[2]);
        final String OUTPUT_FILE = args[3];
        
        final int N_COORDINATE_COLUMNS = 2;
        final int N_COLUMNS = N_COORDINATE_COLUMNS+N_SAMPLES_TOTAL*BAMtracks.N_SIGNALS;
        
        int i, j;
        int length, nSamples;
        double delS, delE, ins, clipL, clipR;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
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
                delS=Double.parseDouble(tokens[i+6]);
                delE=Double.parseDouble(tokens[i+7]);
                ins=Double.parseDouble(tokens[i+8]);
                clipL=Double.parseDouble(tokens[i+9]);
                clipR=Double.parseDouble(tokens[i+10]);
                if (delS!=0 || delE!=0 || ins!=0 || clipL!=0 || clipR!=0) {
                    for (j=0; j<BAMtracks.N_SIGNALS; j++) bw.write(tokens[i+j]+",");
                    bw.newLine();
                }
                nSamples++;
                if (nSamples==N_SAMPLES_OUTPUT) break;
            }
            str=br.readLine();
        }
        br.close(); bw.close();
    }

}