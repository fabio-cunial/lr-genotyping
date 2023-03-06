import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;


/**
 * chr21 length=46709983
 * chr22 length=50818468
 */
public class Pigv {
    /**
     * Internal CIGAR codes
     */
    private static final byte CIGAR_EMPTY = 0;
    private static final byte CIGAR_MATCH = 1;
    private static final byte CIGAR_MISMATCH = 2;
    private static final byte CIGAR_DELETION = 3;
    private static final byte CIGAR_INSERTION = 4;
    
    /**
     * Colormap
     */
    private static final int COLOR_EMPTY = 0x0006080C;
    private static final int COLOR_MATCH = 0x00121415;
    private static final int COLOR_MISMATCH = 0x00829192;
    private static final int COLOR_DELETION = 0x00447194;
    private static final int COLOR_INSERTION = 0x0092182D;
    private static final int COLOR_BACKGROUND = COLOR_EMPTY;
    
    /**
     * Samples divisor
     */
    private static final int DIVISOR_HEIGHT_PIXELS = 1;
    private static final int COLOR_DIVISOR = 0x00DFE1E2;
    
    /**
     * Default sizes of array $cigars$.
     */
    private static final int COVERAGE_PER_SAMPLE = 10;
    private static final int PIXELS_PER_POSITION = 3;
    
    /**
     * Properties of the SV
     */
    private static int svStart, svLength;
    private static String chromosome, pass, svType;
    
    /**
     * Properties of the aligments
     */
    private static int frameStart, frameEnd;
    private static Sample[] samples;
    private static int nSamples;
    
    
    
    /**
     * @param args 0: VCF file with just the name of each column and a single
     * SV call.
     */
    public static void main(String[] args) throws IOException {
        final String VCF_FILE = args[0];
        final String SAM_DIR = args[1];
        final int HORIZONTAL_SLACK = Integer.parseInt(args[2]);
        final int CHROMOSOME_LENGTH = Integer.parseInt(args[3]);
        final String OUTPUT_FILE = args[4];
        
        int i, j, k, x, y;
        int fromX, fromY, imageWidth, imageHeight;
        String str, header, sv, pass;
        BufferedReader br;
        BufferedImage image;
        String[] tokens_header, tokens_sv;
        
        System.err.println("Loading VCF file...");
        br = new BufferedReader(new FileReader(VCF_FILE));
        header=br.readLine(); sv=br.readLine();
        br.close();
        tokens_header=header.split("\t");
        nSamples=tokens_header.length-9;
        samples = new Sample[nSamples];
        tokens_sv=sv.split("\t");
        frameStart=Integer.parseInt(tokens_sv[1])-1-HORIZONTAL_SLACK;
        if (frameStart<0) frameStart=0;
        chromosome=tokens_sv[0];
        svStart=Integer.parseInt(tokens_sv[1])-1;
        pass=tokens_sv[6];
        svType=getField(tokens_sv[7],SVTYPE_STR);
        frameEnd=-1; svLength=-1;
        str=getField(tokens_sv[7],SVLEN_STR);
        if (str!=null) {
            svLength=Integer.parseInt(str);
            frameEnd=svStart+svLength+HORIZONTAL_SLACK;
        }
        if (frameEnd==-1) {
            str=getField(tokens_sv[7],END_STR);
            if (str!=null) frameEnd=Integer.parseInt(str)-1+HORIZONTAL_SLACK;
        }
        if (frameEnd==-1) frameEnd=svStart+HORIZONTAL_SLACK;
        if (frameEnd>=CHROMOSOME_LENGTH) frameEnd=CHROMOSOME_LENGTH-1;    
        
        System.err.print("Loading SAM files... ");
        for (i=0; i<nSamples; i++) {
            samples[i] = new Sample();
            samples[i].id=tokens_header[9+i];
            samples[i].genotype=tokens_sv[9+i].substring(0,3);
            samples[i].loadCigars(SAM_DIR+"/"+samples[i].id+".sam",frameStart,frameEnd);
        }
        System.err.println(" done. "+nSamples+" samples loaded.");
        
        System.err.println("Drawing image...");
        imageHeight=0;
        for (i=0; i<nSamples; i++) imageHeight+=samples[i].cigars_last+1+DIVISOR_HEIGHT_PIXELS;
        imageWidth=(frameEnd-frameStart+1)*PIXELS_PER_POSITION;
        image = new BufferedImage(imageWidth,imageHeight,BufferedImage.TYPE_INT_RGB);
        System.err.println("imageWidth="+imageWidth+" imageHeight="+imageHeight);
		for (x=0; x<imageWidth; x++) {
		    for (y=0; y<imageHeight; y++) image.setRGB(x,y,COLOR_BACKGROUND);
		}
        Sample.order=Sample.ORDER_GENOTYPE;
        Arrays.sort(samples);
        fromX=0; fromY=0;
        for (i=0; i<nSamples; i++) {
            fromY=samples[i].draw(fromX,fromY,image,imageWidth,imageHeight)+1;
            for (j=0; j<DIVISOR_HEIGHT_PIXELS; j++) {
                for (k=0; k<imageWidth; k++) image.setRGB(k,fromY,COLOR_DIVISOR);
                fromY++;
            }
        }
        ImageIO.write(image,"png",new File(OUTPUT_FILE));
    }
    
    
    private static final int nHaplotypes(String genotype) {
        return (genotype.charAt(0)==1?1:0) + (genotype.charAt(2)==1?1:0);
    }
    
    
    /**
     * An individual at a variant
     */
    public static class Sample implements Comparable {
        /**
         * Global sorting criterion. The sample that is most different from the
         * reference always comes first.
         */
        public static final byte ORDER_GENOTYPE = 0;
        public static final byte ORDER_CIGARS = 1;
        public static byte order;
        
        /**
         * Information about the individual
         */
        public String id, genotype;
        
        /**
         * Dimensions: 0=read, 1=position in the current frame.
         */
        public byte[][] cigars;
        public int cigars_nColumns, cigars_last;
        
        
        public Sample() { }
        
        
        /**
         * @param samFile assumed to contain no header;
         * @param frameStart first position in the current frame (zero-based).
         */
        public final void loadCigars(String samFile, int frameStart, int frameEnd) throws IOException {
            boolean isEmpty;
            int i, j, k;
            String str;
            BufferedReader br;
            String[] tokens;
            
            // Loading all cigar strings
            cigars_nColumns=(frameEnd-frameStart+1)*PIXELS_PER_POSITION;
            cigars = new byte[COVERAGE_PER_SAMPLE][cigars_nColumns];
            cigars_last=-1;
            br = new BufferedReader(new FileReader(samFile));
            str=br.readLine();
            while (str!=null) {
                tokens=str.split("\t");
                loadCigar(tokens[5],Integer.parseInt(tokens[3])-1,frameStart);
                str=br.readLine();
            }
            br.close();
            
            // Removing empty rows
            k=-1;
            for (i=0; i<=cigars_last; i++) {
                isEmpty=true;
                for (j=0; j<cigars_nColumns; j++) {
                    if (cigars[i][j]!=CIGAR_EMPTY) {
                        isEmpty=false;
                        break;
                    }
                }
                if (!isEmpty) {
                    k++;
                    if (k!=i) System.arraycopy(cigars[i],0,cigars[k],0,cigars_nColumns);
                }
            }
            cigars_last=k;
        }
        
        
        /**
         * Appends a new row to matrix $cigars$. 
         *
         * @param cigarStart first position of the reference used by $cigar$
         * (zero-based);
         * @param frameStart first position in the current frame (zero-based).
         */
        private final void loadCigar(String cigar, int cigarStart, int frameStart) {
            boolean previousInsertion;
    		char c, d;
    		int i, j, k, p, q;
            int from, length;
    		final int cigarLength = cigar.length();
	        
            // Allocating a new row in $cigars$.
            cigars_last++;
            if (cigars_last==cigars.length) {
                byte[][] newArray = new byte[cigars.length<<1][cigars_nColumns];
                for (i=0; i<cigars.length; i++) System.arraycopy(cigars[i],0,newArray[i],0,cigars_nColumns);
                cigars=newArray;
            }
            Arrays.fill(cigars[cigars_last],CIGAR_EMPTY);
            
            // Skipping the initial hard and soft clips, if any.
    		from=-1;
    		for (i=0; i<cigarLength; i++) {
                c=cigar.charAt(i);
    			if (c=='H') {
                    from=i;
            		for (j=i+1; j<cigarLength; j++) {
                        d=cigar.charAt(j);
            			if (d=='S') {
                            from=j;
                            break;
                        }
                        else if (!Character.isDigit(d)) break;
            		}
                    break;
                }
                else if (c=='S') {
                    from=i;
                    break;
                }
                else if (!Character.isDigit(c)) break;
    		}
            from++;
    		
            // Loading the rest
    		p=from; q=(cigarStart-frameStart)*PIXELS_PER_POSITION; previousInsertion=false;
    		for (i=from; i<cigarLength; i++) {
    			c=cigar.charAt(i);
    			if (c=='M' || c=='=') {
    				length=Integer.parseInt(cigar.substring(p,i));
                    for (j=0; j<length; j++) {
                        if (q>=0) cigars[cigars_last][q]=previousInsertion?CIGAR_INSERTION:CIGAR_MATCH;
                        q++;
                        previousInsertion=false;
                        for (k=1; k<PIXELS_PER_POSITION && q<cigars_nColumns; k++) {
                            if (q>=0) cigars[cigars_last][q]=CIGAR_MATCH;
                            q++;
                        }
                        if (q>=cigars_nColumns) break;
                    }
                    p=i+1;
    			}
    			else if (c=='X') {
    				length=Integer.parseInt(cigar.substring(p,i));
                    for (j=0; j<length; j++) {
                        if (q>=0) cigars[cigars_last][q]=previousInsertion?CIGAR_INSERTION:CIGAR_MISMATCH;
                        q++;
                        previousInsertion=false;
                        for (k=1; k<PIXELS_PER_POSITION && q<cigars_nColumns; k++) {
                            if (q>=0) cigars[cigars_last][q]=CIGAR_MISMATCH;
                            q++;
                        }
                        if (q>=cigars_nColumns) break;
                    }
                    p=i+1;
    			}
    			else if (c=='D' || c=='N') {
    				length=Integer.parseInt(cigar.substring(p,i));
                    for (j=0; j<length; j++) {
                        if (q>=0) cigars[cigars_last][q]=previousInsertion?CIGAR_INSERTION:CIGAR_DELETION;
                        q++;
                        previousInsertion=false;
                        for (k=1; k<PIXELS_PER_POSITION && q<cigars_nColumns; k++) {
                            if (q>=0) cigars[cigars_last][q]=CIGAR_DELETION;
                            q++;
                        }
                        if (q>=cigars_nColumns) break;
                    }
                    p=i+1;
    			}
    			else if (c=='I') {
    				length=Integer.parseInt(cigar.substring(p,i));
                    if (q-1>=0) cigars[cigars_last][q-1]=CIGAR_INSERTION;
                    previousInsertion=true;
                    p=i+1;
    			}
    			else if (c=='S' || c=='H') break;
                if (q>=cigars_nColumns) break;
    		}
        }
        
        
        public int compareTo(Object other) {
            Sample otherSample = (Sample)other;
            
            if (order==ORDER_GENOTYPE) {
                final int n = nHaplotypes(genotype);
                final int otherN = nHaplotypes(otherSample.genotype);
                if (n>otherN) return -1;
                else if (n<otherN) return 1;
                return 0;
            }
            else if (order==ORDER_CIGARS) {
                final double n = mismatchRate();
                final double otherN = otherSample.mismatchRate();
                if (n>otherN) return -1;
                else if (n<otherN) return 1;
                return 0;
            }
            else return 0;
        }
        
        
        private final double mismatchRate() {
            int i, j;
            double numerator, denominator;
            
            numerator=0.0; denominator=0.0;
            for (i=0; i<=cigars_last; i++) {
                for (j=0; j<cigars_nColumns; j++) {
                    if (cigars[i][j]==CIGAR_EMPTY) continue;
                    denominator+=1.0;
                    if (cigars[i][j]!=CIGAR_MATCH) numerator+=1.0;
                }
            }
            return numerator/denominator;
        }
        
        
        /**
         * @param array1,array2 assumed to be of the same length.
         */
        private static final double similarity(byte[] array1, byte[] array2) {
            int i;
            final int length = array1.length;
            double numerator;
            
            numerator=0.0;
            for (i=0; i<length; i++) {
                if (array1[i]==array2[i]) numerator+=1.0;
            }
            return numerator/length;
        }
        
        
        /**
         * @param imageWidth,imageHeight in pixels;
         * @return the last Y coordinate used by the drawing.
         */
        private final int draw(int fromX, int fromY, BufferedImage image, int imageWidth, int imageHeight) {
            int x, y;
            int color;
            
            for (y=0; y<=cigars_last; y++) {
                if (fromY+y>=imageHeight) break;
                for (x=0; x<cigars_nColumns; x++) {                               
                    if (fromX+x>=imageWidth) break;
                    color=COLOR_EMPTY;
                    if (cigars[y][x]==CIGAR_MATCH) color=COLOR_MATCH;
                    else if (cigars[y][x]==CIGAR_MISMATCH) color=COLOR_MISMATCH;
                    else if (cigars[y][x]==CIGAR_DELETION) color=COLOR_DELETION;
                    else if (cigars[y][x]==CIGAR_INSERTION) color=COLOR_INSERTION;
                    image.setRGB(fromX+x,fromY+y,color);
                }
            }
            return fromY+cigars_last;
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
    
    
	/**
	 * @return -1 iff the type cannot be determined.
	 */
	private static final int getType_infoField(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return TYPE_DELETION;
		else if (type.equalsIgnoreCase(DEL_INV_STR)) return TYPE_DEL_INV;
		else if ( type.equalsIgnoreCase(INS_STR) || 
			      type.equalsIgnoreCase(INS_ME_STR) || 
				  type.equalsIgnoreCase(INS_NOVEL_STR)
				) return TYPE_INSERTION;
		else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return TYPE_DUPLICATION;
		else if (type.equalsIgnoreCase(INV_STR)) return TYPE_INVERSION;
		else if (type.equalsIgnoreCase(INV_DUP_STR)) return TYPE_INV_DUP;
		else if (type.equalsIgnoreCase(CNV_STR)) return TYPE_CNV;
		else if (type.equalsIgnoreCase(BND_STR)) return TYPE_BREAKEND;
		else if (type.equalsIgnoreCase(TRA_STR)) return TYPE_TRANSLOCATION;
		else return -1;
	}

}