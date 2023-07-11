import java.util.PriorityQueue;
import java.io.*;

/**
 * Writes to a VCF all the INS and DEL that are in HW equilibrium, and to
 * another VCF all the INS and DEL that are not in HW equilibrium.
 *
 * Remark: MIN_CHI_SQUARE=16.5 has p-value=0.000049.
 */
public class HardyWeinbergSVs {

	public static void main(String[] args) throws IOException {
	    final String VCF_FILE = args[0];
		final double MIN_CHI_SQUARE = Double.parseDouble(args[1]);
        final int TOP_K = Integer.parseInt(args[2]);
        
        int i, n, a, b;
        int nSamples, length, n11, n01, n00;
        double p, q, e11, e01, e00, chiSquare;
        String str, info;
        BufferedReader br;
        BufferedWriter inHW, notinHW;
        SV tmpSV;
        PriorityQueue<SV> queue;
        String[] tokens;
        
        queue = new PriorityQueue<SV>();
        inHW = new BufferedWriter(new FileWriter(VCF_FILE+".inHW.vcf"));
        notinHW = new BufferedWriter(new FileWriter(VCF_FILE+".notinHW.vcf"));
        br = new BufferedReader(new FileReader(VCF_FILE));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                inHW.write(str+"\n");
                notinHW.write(str+"\n");
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            info=tokens[7].toLowerCase();
            if (info.indexOf("svtype=ins")<0 && info.indexOf("svtype=del")<0) {
                str=br.readLine();
                continue;
            }
            n=0; n11=0; n01=0; n00=0;
            for (i=9; i<tokens.length; i++) {
                if (tokens[i].indexOf("0/0")>=0) { n++; n00++; }
                else if (tokens[i].indexOf("0/1")>=0 || tokens[i].indexOf("1/0")>=0) { n++; n01++; }
                else if (tokens[i].indexOf("1/1")>=0) { n++; n11++; }
            }
            if (n<tokens.length-9) System.err.println("Some samples miss a genotype:  "+str);
            nSamples=n11+n01+n00;
            if (nSamples==0) {
                System.err.println("No genotype at all:  "+str);
                str=br.readLine();
                continue;
            }
            p=(2.0*n11+n01)/(2.0*nSamples); q=1-p;
            e11=p*p*nSamples; e01=2*p*q*nSamples; e00=q*q*nSamples;
            chiSquare=((n11-e11)*(n11-e11))/e11 + ((n01-e01)*(n01-e01))/e01 + ((n00-e00)*(n00-e00))/e00;
            if (chiSquare>=MIN_CHI_SQUARE) {
                notinHW.write(str+"\n");
                queue.add(new SV(str,chiSquare));
            }
            else inHW.write(str+"\n");
            str=br.readLine();
        }
        br.close(); inHW.close(); notinHW.close();
        
        // Printing $TOP_K$ SVs that are not in HW.
        n=TOP_K>queue.size()?queue.size():TOP_K;
        for (i=0; i<n; i++) {
            tmpSV=queue.poll();
            tokens=tmpSV.str.split("\t");
            info=tokens[7].toLowerCase();
            if (info.indexOf("svtype=ins")>=0) System.out.println(tmpSV.chiSquare+"\t"+tokens[0]+"\t"+tokens[1]+"\t"+tokens[1]+"\t"+tmpSV.str);
            else if (info.indexOf("svtype=del")>=0) {
                a=tokens[7].indexOf("SVLEN="); b=tokens[7].indexOf(";",a+6);
                if (b<0) b=tokens[7].length();
                length=Integer.parseInt(tokens[7].substring(a+6,b));
                if (length<0) length=-length;
                System.out.println(tmpSV.chiSquare+"\t"+tokens[0]+"\t"+tokens[1]+"\t"+(Integer.parseInt(tokens[1])+length-1)+"\t"+tmpSV.str);
            }
        }
	}
    
    
    private static class SV implements Comparable {
        public String str;
        public double chiSquare;
        
        public SV(String s, double c) {
            this.str=s; this.chiSquare=c;
        }
        
        public int compareTo(Object other) {
            SV otherSV = (SV)other;
            if (chiSquare<otherSV.chiSquare) return 1;
            else if (chiSquare>otherSV.chiSquare) return -1;
            return 0;
        }
    }
    
}
