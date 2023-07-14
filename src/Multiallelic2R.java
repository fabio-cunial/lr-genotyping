import java.util.HashMap;
import java.io.*;

/**
 * 
 *
 * The R script can then be executed by doing:
 *
 * Rscript out.r > stats.txt
 */
public class Multiallelic2R {

	public static void main(String[] args) throws IOException {
	    final String VCF_FILE = args[0];
        final int N_PERMUATIONS = Integer.parseInt(args[1]);
        final String R_FILE = args[2];
        
        int i;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        Genotype query;
        Genotype tmpGenotype = new Genotype();
        HashMap<Genotype,Genotype> map;
        String[] tokens;
        Genotype[] tmpArray;
        
        map = new HashMap<Genotype,Genotype>();
        br = new BufferedReader(new FileReader(VCF_FILE));
        bw = new BufferedWriter(new FileWriter(R_FILE));
        bw.write("library(HardyWeinberg)\n");
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t"); map.clear();
            for (i=9; i<tokens.length; i++) {
                tmpGenotype.load(tokens[i]);
                query=map.get(tmpGenotype);
                if (query!=null) query.count++;
                else {
                    query=tmpGenotype.clone();
                    query.count=1;
                    map.put(query,query);
                }
            }
            tmpArray = new Genotype[map.size()];
            map.keySet().toArray(tmpArray);
            toRmatrix(tmpArray,N_PERMUATIONS,bw);
            str=br.readLine();
        }
        br.close(); bw.close();
	}
    
    
    private static final void toRmatrix(Genotype[] array, int nPermutations, BufferedWriter bw) throws IOException {
        int i, j;
        int max, nRows;
        final int length = array.length;
        int[][] matrix;
        
        max=0;
        for (i=0; i<length; i++) {
            if (array[i].hap1>max) max=array[i].hap1;
        }
        nRows=max+1;
        matrix = new int[nRows][nRows];
        for (i=0; i<length; i++) matrix[array[i].hap1][array[i].hap2]=array[i].count;
        bw.write("A = matrix( \n c(");
        for (i=0; i<nRows-1; i++) {
            for (j=0; j<nRows; j++) bw.write(matrix[i][j]+",");
        }
        bw.write(matrix[nRows-1][0]+"");
        for (j=1; j<nRows; j++) bw.write(","+matrix[nRows-1][j]);
        bw.write("), \n nrow = "+nRows+",\n ncol = "+nRows+",\n byrow = TRUE\n )\n");
        bw.write("rownames(A) = c(");
        bw.write("\""+0+"\"");
        for (i=1; i<nRows; i++) bw.write(",\""+i+"\"");
        bw.write(") \n colnames(A) = c(");
        bw.write("\""+0+"\"");
        for (i=1; i<nRows; i++) bw.write(",\""+i+"\"");
        bw.write(") \n");
        bw.write("out <- HWPerm.mult(A,nperm="+nPermutations+")\n");
    }
    
    
    private static class Genotype implements Comparable {
        public int hap1, hap2;  // $hap1>=hap2$
        public int count;
        
        public Genotype() { hap1=-1; hap2=-1; count=-1; }
        
        public Genotype(int h1, int h2) {
            hap1=h1; hap2=h2; count=0;
        }
        
        public void load(String str) {
            int p, q, h1, h2;
            String hap;
            
            p=str.indexOf("/");
            if (p<0) p=str.indexOf("|");
            q=str.indexOf(":",p+1);
            hap=str.substring(0,p);
            h1=hap.equals(".")?0:Integer.parseInt(hap);
            hap=str.substring(p+1,q);
            h2=hap.equals(".")?0:Integer.parseInt(hap);
            if (h1>h2) { hap1=h1; hap2=h2; }
            else { hap1=h2; hap2=h1; }
        }
        
        public Genotype clone() {
            return new Genotype(hap1,hap2);
        }
        
        public int compareTo(Object other) {
            Genotype otherGenotype = (Genotype)other;
            if (hap1<otherGenotype.hap1) return 1;
            else if (hap1>otherGenotype.hap1) return -1;
            if (hap2<otherGenotype.hap2) return 1;
            else if (hap2>otherGenotype.hap2) return -1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Genotype otherGenotype = (Genotype)other;
            return hap1==otherGenotype.hap1 && hap2==otherGenotype.hap2;
        }
        
        public int hashCode() {
            return (hap1+""+hap2).hashCode();
        }
    }
    
}
