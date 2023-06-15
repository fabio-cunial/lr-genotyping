import java.io.*;


/**
 * Prints histograms on the following sniffles2 tags:
 *
 * SUPPORT            Number of reads supporting the structural variation.
 * SUPPORT_INLINE     Number of reads supporting an INS/DEL SV (non-split events 
 *                    only).
 * SUPPORT_LONG       Number of soft-clipped reads putatively supporting the
 *                    long insertion.
 * CONSENSUS_SUPPORT  Number of reads that support the generated insertion (INS)
 *                    consensus sequence.
 * STDEV_POS          Standard deviation of structural variation start position.
 * STDEV_LEN          Standard deviation of structural variation length.
 */
public class Sniffles2Histograms {
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String VCF_FILE = args[0];
        
        int i, j;
        String str;
        BufferedReader br;
        int[][] histograms;
        String[] tokens;
        
        histograms = new int[31][6];
        br = new BufferedReader(new FileReader(VCF_FILE));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            str=getField(tokens[7],"SUPPORT");
            if (str!=null) {
                i=Integer.parseInt(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][0]++;
            }
            str=getField(tokens[7],"SUPPORT_INLINE");
            if (str!=null) {
                i=Integer.parseInt(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][1]++;
            }
            str=getField(tokens[7],"SUPPORT_LONG");
            if (str!=null) {
                i=Integer.parseInt(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][2]++;
            }
            str=getField(tokens[7],"CONSENSUS_SUPPORT");
            if (str!=null) {
                i=Integer.parseInt(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][3]++;
            }
            str=getField(tokens[7],"STDEV_POS");
            if (str!=null) {
                i=(int)Double.parseDouble(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][4]++;
            }
            str=getField(tokens[7],"STDEV_LEN");
            if (str!=null) {
                i=(int)Double.parseDouble(str);
                if (i>=histograms.length) i=histograms.length-1;
                histograms[i][5]++;
            }
            str=br.readLine();
        }
        br.close();
        for (i=0; i<histograms.length; i++) {
            for (j=0; j<histograms[i].length; j++) System.out.print(histograms[i][j]+",");
            System.out.println();
        }
    }
    
    
	/**
	 * Basic constants
	 */
	public static final char COMMENT = '#';
	public static final String SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String PRECISE_STR = "PRECISE";
	public static final String IMPRECISE_STR = "IMPRECISE";
	public static final String END_STR = "END";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
	public static final String CHR2_STR = "CHR2";
	public static final String CT_STR = "CT";
	public static final String CT_325_STR = "'3to5'";
	public static final String CT_523_STR = "'5to3'";
	public static final String CT_525_STR = "'5to5'";
	public static final String CT_323_STR = "'3to3'";
    
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
	public static final String INS_NOVEL_STR = "INS:NOVEL";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";
    
	/**
     *
	 */
	public static final byte TYPE_INSERTION = 1;
	public static final byte TYPE_DELETION = 2;
	public static final byte TYPE_DEL_INV = 3;
	public static final byte TYPE_INVERSION = 4;
	public static final byte TYPE_INV_DUP = 5;
	public static final byte TYPE_DUPLICATION = 6;
	public static final byte TYPE_CNV = 7;
	public static final byte TYPE_BREAKEND = 8;
	public static final byte TYPE_TRANSLOCATION = 9;
    
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which our programs disregard for 
	 * simplicity. Some callers use CILEN to express a "confidence interval 
	 * around inserted/deleted material between breakends": we interpret CILEN
	 * exactly like CIEND, and we ignore it for insertions (since representing
	 * variable-length insertion strings complicates our code).
	 */
	public static final String CI_SEPARATOR = ",";
	public static final String CIPOS_STR = "CIPOS";
	public static final String CIEND_STR = "CIEND";
	public static final String STD_START1_STR = "STD_quant_start";
	public static final String STD_START2_STR = "STD_POS1";
	public static final String STD_END1_STR = "STD_quant_stop";
	public static final String STD_END2_STR = "STD_POS2";
	public static final String CILEN_STR = "CILEN";
	public static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
	public static final int CIEND_STR_LENGTH = CIEND_STR.length();
	public static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
	public static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
	public static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
	public static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
	public static final int CILEN_STR_LENGTH = CILEN_STR.length();
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	public static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}
    
}