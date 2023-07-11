import java.util.Arrays;
import java.io.*;


/**
 * chr21 length: 46709983
 * chr22 length: 50818468
 */
public class Histogram {
    
	public static void main(String[] args) throws IOException {
	    final String VCF_IN_HW = args[0];
        final String VCF_NOTIN_HW = args[1];
        final String CHR = args[2];
        final int CHR_LENGTH = Integer.parseInt(args[3]);
		final int BIN_LENGTH_POSITION = Integer.parseInt(args[4]);
        final int BIN_LENGTH_LENGTH = Integer.parseInt(args[5]);
        
        int i, j;
        final int nRows_position = CHR_LENGTH/BIN_LENGTH_POSITION+1;
        final int nRows_length = 10000/BIN_LENGTH_LENGTH+1;
        final int nColumns = 4;
        BufferedReader br;
        int[][] histogram_position, histogram_length;
        
        histogram_position = new int[nRows_position][nColumns];
        for (i=0; i<nRows_position; i++) Arrays.fill(histogram_position[i],0);
        histogram_length = new int[nRows_length][nColumns];
        for (i=0; i<nRows_length; i++) Arrays.fill(histogram_length[i],0);
        br = new BufferedReader(new FileReader(VCF_IN_HW));
        buildHistograms(histogram_position,BIN_LENGTH_POSITION,histogram_length,BIN_LENGTH_LENGTH,0,br,CHR);
        br.close();
        br = new BufferedReader(new FileReader(VCF_NOTIN_HW));
        buildHistograms(histogram_position,BIN_LENGTH_POSITION,histogram_length,BIN_LENGTH_LENGTH,2,br,CHR);
        br.close();
        for (i=0; i<nRows_position; i++) {
            for (j=0; j<nColumns; j++) System.out.print(histogram_position[i][j]+",");
            System.out.println();
        }
        for (i=0; i<nRows_length; i++) {
            for (j=0; j<nColumns; j++) System.err.print(histogram_length[i][j]+",");
            System.err.println();
        }
	}
    
    
    /**
     * @param histogram_position a table whose rows are position bins of length
     * $binLength_position$ each, and whose columns are:
     * nINS in HW equilibrium
     * nDEL in HW equilibrium
     * nINS not in HW equilibrium
     * nDEL not in HW equilibrium;
     * @param histogram_length a table with the same columns as above, but whose
     * rows are SV length bins of length $binLength_length$ each.
     */
    private static final void buildHistograms(int[][] histogram_position, int binLength_position, int[][] histogram_length, int binLength_length, int firstColumn, BufferedReader br, String chr) throws IOException {
        int p, q;
        int bin_position, bin_length;
        String str, info;
        String[] tokens;
        
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            if (!tokens[0].equalsIgnoreCase(chr)) {
                str=br.readLine();
                continue;
            }
            bin_position=(Integer.parseInt(tokens[1])-1)/binLength_position;
            info=tokens[7].toLowerCase();
            p=info.indexOf("svlen="); q=info.indexOf(";",p+6);
            if (q<0) q=info.length();
            bin_length=Integer.parseInt(info.substring(p+6,q))/binLength_length;
            if (bin_length<0) bin_length=-bin_length;
            if (bin_length>histogram_length.length-1) bin_length=histogram_length.length-1;
            if (info.indexOf("svtype=ins")>=0) {
                histogram_position[bin_position][firstColumn]++;
                histogram_length[bin_length][firstColumn]++;
            }
            else if (info.indexOf("svtype=del")>=0) {
                histogram_position[bin_position][firstColumn+1]++;
                histogram_length[bin_length][firstColumn+1]++;
            }
            str=br.readLine();
        }
        br.close();
    }    
    
}
