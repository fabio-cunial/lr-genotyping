import java.io.*;
import java.awt.*;
import java.awt.image.*;
import javax.imageio.*;


/** 
 * Draws an image with a row for every sample. Prints a histogram with the sum
 * of values over all samples.
 *
 * Remark: before visualizing the file, it must be cleaned as follows:
 * sed 's/,,/,/g' out.txt > out-cleaned.txt
 *
 * Remark: after visualization, on can concatenate all images as follows (imagemagick):
 * convert coverage.png clip.png del.png ins.png kmer.png -append all.png
 */
public class BAMtracksPrint {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int SIGNAL = Integer.parseInt(args[1]);
        final boolean USE_LOG = Integer.parseInt(args[2])==1;
        final String OUTPUT_FILE_IMAGE = args[3];
        final String OUTPUT_FILE_HISTOGRAM = args[4];
        
        int i, x, y, n;
        int nRows, nColumns, length, maxRow, maxColumn;
        double value, max, sum, average, stddev;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        BufferedImage image;
        String[] tokens;
        
        System.err.println("Allocating image...");
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        tokens=str.split(",");
        nRows=(tokens.length-2)/BAMtracks.N_SIGNALS;
        nColumns=0; max=0.0; maxRow=-1; maxColumn=-1;
        while (str!=null) {
            nColumns++;
            tokens=str.split(",");
            length=tokens.length;
            for (i=2+SIGNAL; i<length; i+=BAMtracks.N_SIGNALS) {
                if (tokens[i].length()==0) continue;
                value=Double.parseDouble(tokens[i]);
                if (value>max) { max=value; maxRow=(i-2)/BAMtracks.N_SIGNALS; maxColumn=nColumns-1; }
            }
            str=br.readLine();
        }
        br.close();
        image = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
        System.err.println("nRows="+nRows+" nColumns="+nColumns+" max="+max+" maxRow="+maxRow+" maxColumn="+maxColumn);
        
        System.err.println("Drawing image...");
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_HISTOGRAM));
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        x=-1;
        while (str!=null) {
            x++;
            tokens=str.split(","); length=tokens.length; sum=0.0; n=0;
            for (i=2+SIGNAL; i<length; i+=BAMtracks.N_SIGNALS) {
                n++;
                if (tokens[i].length()==0) {
                    image.setRGB(x,(i-2)/BAMtracks.N_SIGNALS,COLOR_EMPTY_VALUE);
                    continue;
                }
                value=Double.parseDouble(tokens[i]);
                sum+=value;
                image.setRGB(x,(i-2)/BAMtracks.N_SIGNALS,getColor(value,max,USE_LOG));
            }
            average=sum/n; stddev=0.0;
            for (i=2+SIGNAL; i<length; i+=BAMtracks.N_SIGNALS) {
                if (tokens[i].length()==0) stddev+=average*average;
                else {
                    value=Double.parseDouble(tokens[i]);
                    stddev+=(value-average)*(value-average);
                }
            }
            bw.write(sum+","+average+","+Math.sqrt(stddev)); bw.newLine();
            str=br.readLine();
        }
        br.close(); bw.close();
        ImageIO.write(image,"png",new File(OUTPUT_FILE_IMAGE));
    }
    
    
    private static final int COLOR_EMPTY_VALUE = 0x0059A67D;
    private static final int COLOR_FLAG_1 = 0x00E61A1A;
    private static final int COLOR_FLAG_2 = 0x00267BD9;
    
    
    private static final int getColor(double value, double max, boolean useLog) {
        int p, q;
        final int MASK = 0x000000FF;
        
        if (value==-2) return COLOR_FLAG_2;
        else if (value==-1) return COLOR_FLAG_1;
        else {
            if (useLog) p=(int)(255*(Math.log(value)/Math.log(max)));
            else p=(int)((255*value)/max);
            p&=MASK;
            q=0x00000000;
            q|=p; q<<=8;
            q|=p; q<<=8;
            q|=p;
            return q;
        }
    }

}