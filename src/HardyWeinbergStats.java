import java.io.*;

/**
 * Reads from STDIN
 */
public class HardyWeinbergStats {

	public static void main(String[] args) throws IOException {
        boolean isIns, isDel;
        int i, n, a, b;
        int nSamples, length, n11, n01, n00;
        double p, q, e11, e01, e00, chiSquare;
        String str, info;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(System.in));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            info=tokens[7].toLowerCase();
            isIns=info.indexOf("svtype=ins")>=0;
            isDel=info.indexOf("svtype=del")>=0;
            if (!isIns && !isDel) {
                str=br.readLine();
                continue;
            }
            n=0; n11=0; n01=0; n00=0;
            for (i=9; i<tokens.length; i++) {
                if (tokens[i].indexOf("0/0")>=0) { n++; n00++; }
                else if (tokens[i].indexOf("0/1")>=0 || tokens[i].indexOf("1/0")>=0) { n++; n01++; }
                else if (tokens[i].indexOf("1/1")>=0) { n++; n11++; }
            }
            if (n<tokens.length-9) System.out.println("Some samples miss a genotype:  "+str);
            nSamples=n11+n01+n00;
            if (nSamples==0) {
                System.out.println("No genotype at all:  "+str);
                str=br.readLine();
                continue;
            }
            p=(2.0*n11+n01)/(2.0*nSamples); q=1-p;
            e11=p*p*nSamples; e01=2*p*q*nSamples; e00=q*q*nSamples;
            chiSquare=((n11-e11)*(n11-e11))/e11 + ((n01-e01)*(n01-e01))/e01 + ((n00-e00)*(n00-e00))/e00;
            a=info.indexOf("svlen="); b=info.indexOf(";",a+6);
            if (b<0) b=info.length();
            length=Integer.parseInt(info.substring(a+6,b));
            if (length<0) length=-length;
            System.out.println(tokens[0]+":"+tokens[1]+" "+(isIns?"INS":"DEL")+" SVLEN="+length+" 0/0="+n00+":"+(int)e00+" 0/1="+n01+":"+(int)e01+" 1/1="+n11+":"+(int)e11+" chiSquare="+chiSquare+"    e00="+e00+" e01="+e01+" e11="+e11+"          ((n11-e11)*(n11-e11))/e11="+(((n11-e11)*(n11-e11))/e11)+"      ((n01-e01)*(n01-e01))/e01="+(((n01-e01)*(n01-e01))/e01)+"       ((n00-e00)*(n00-e00))/e00="+(((n00-e00)*(n00-e00))/e00));
            str=br.readLine();
        }
        br.close();
	}
    
}
