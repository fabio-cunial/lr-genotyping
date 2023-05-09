import java.io.*;
import java.awt.*;
import java.awt.image.*;
import javax.imageio.*;


/**
 * 
 */
public class BAMtracksPrint {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int SIGNAL = Integer.parseInt(args[1]);
        final String OUTPUT_FILE = args[2];
        
        int i, x, y;
        int nRows, nColumns, length;
        double value, max;
        String str;
        BufferedReader br;
        BufferedImage image;
        String[] tokens;
        
        System.err.println("Allocating image...");
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        tokens=str.split(",");
        nRows=(tokens.length-2)/5;
        nColumns=0; max=0.0;
        while (str!=null) {
            nColumns++;
            tokens=str.split(",");
            length=tokens.length;
            for (i=2+SIGNAL; i<length; i+=5) {
                if (tokens[i].length()==0) continue;
                value=Double.parseDouble(tokens[i]);
                if (value>max) max=value;
            }
            str=br.readLine();
        }
        br.close();
        image = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
        System.err.println("nRows="+nRows+" nColumns="+nColumns+" max="+max);
        
        System.err.println("Drawing image...");
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        x=-1;
        while (str!=null) {
            x++;
            tokens=str.split(","); length=tokens.length;
            for (i=2+SIGNAL; i<length; i+=5) {
                if (tokens[i].length()==0) {
                    image.setRGB(x,(i-2)/5,COLOR_EMPTY_VALUE);
                    continue;
                }
                value=Double.parseDouble(tokens[i]);
                image.setRGB(x,(i-2)/5,getColor(value,max));
            }
            str=br.readLine();
        }
        br.close();
        ImageIO.write(image,"png",new File(OUTPUT_FILE));
    }
    
    
    private static final int COLOR_EMPTY_VALUE = 0x0059A67D;
    private static final int COLOR_FLAG_1 = 0x00E61A1A;
    private static final int COLOR_FLAG_2 = 0x00267BD9;
    
    
    private static final int getColor(double value, double max) {
        int p, q;
        final int MASK = 0x000000FF;
        
        if (value==-2) return COLOR_FLAG_2;
        else if (value==-1) return COLOR_FLAG_1;
        else {
            p=(int)(255*(Math.log(value)/Math.log(max)));
            p&=MASK;
            q=0x00000000;
            q|=p; q<<=8;
            q|=p; q<<=8;
            q|=p;
            return q;
        }
    }

}