/*
 * Author: Steve Monger 2020
 * 
 * Description: Merges all variant and annotation information to determine motif changes and exon-intron phase.
 * These are integrated to determine logistic regression probabilities.
 * 
 * Input: GTF annotation
 * Output: One line per gene, chr, geneStart, geneEnd, geneName, strand, donorStartPos(exon1..n-1), DonorEndPos(exon1..n-1), AcceptorStartPos(exon2..n), AcceptorEndPos(exon2..n), donorExonID, acceptorExonID.
 */

//Description: From GTF annotation, output splice site intervals for multi exon genes.
//Exon start/end positions which are also transcript start/end positions are ignored, as are duplicate exons.
//Output order: chr, geneStart, geneEnd, geneName, strand, donorStartPos(exon1..n-1), DonorEndPos(exon1..n-1),
//AcceptorStartPos(exon2..n), AcceptorEndPos(exon2..n), donorExonID, acceptorExonID.
import java.io.*;
import java.util.Arrays;

class getSpliceSiteIntervalsFromGTF
{
    public static int[] accStartPos = new int[10000];
    public static int[] accEndPos = new int[10000];
    public static int[] donStartPos = new int[10000];
    public static int[] donEndPos = new int[10000];
    public static int donIndex = 0;
    public static int accIndex = 0;
    public static String[] donNames = new String[10000];
    public static String[] accNames = new String[10000];
    public static String geneName = "";
    public static int geneStart = Integer.MAX_VALUE;
    public static int geneEnd = 0;
    public static int codingStart = 0;
    public static int codingEnd = 0;
    public static String chr = "";
    public static String strand = "";
    public static String prevLine = "firstLine";
    public static boolean duplicateExon = false;
    public static String[] outPos = new String[50000];
    public static String[] outNeg = new String[50000];
    public static int outPosIndex = 0;
    public static int outNegIndex = 0;

    public static void main(String[] args)
    {	
        String basenameGTF=args[0];
        String line = "";
        try
        {
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            while ((line = in.readLine()) != null)
            {
                if (!line.startsWith("#"))
                {
	                getIntervals(line);
                }
            }
        in.close();
        }
        catch (Exception e)
        {
            System.err.println("Error: ".concat(e.getMessage()));
        }

        writeToFile(basenameGTF);
    }

	public static void getIntervals(String line)
    {
		String[] sep = line.split("\t");
	    if (sep[2].equals("gene"))
        {
        	//remove the final exon end pos from donor output for the previous gene
	        if (prevLine.equals("exon"))
            {
                if (donIndex!=0)
                {
	                donStartPos[donIndex-1]=0;
                    donEndPos[donIndex-1]=0;
                    donNames[donIndex-1]="";
                    donIndex--;
	            }	            
	        }

	        //output previous
	        if (!prevLine.equals("firstLine"))
            {
	            formatOutputAndReset();
	        }

	        //update variables
	        geneName = sep[8].substring(sep[8].indexOf("gene_name \""));
	        geneName = geneName.substring(11, geneName.indexOf("\";"));
	        chr = sep[0];
		    geneStart = Integer.parseInt(sep[3]);
		    geneEnd = Integer.parseInt(sep[4]);
	        strand = sep[6];      
		    prevLine = "gene";
	    }

	    //exclude transcript end positions from splice site intervals
        else if (sep[2].equals("transcript"))
        {
	        if (prevLine.equals("exon"))
            {
                if (donIndex!=0)
                {
	                donStartPos[donIndex-1]=0;
                    donEndPos[donIndex-1]=0;
                    donNames[donIndex-1]="";
                    donIndex--;
	            }	            
	        }
		    prevLine = "transcript";
	    }

        //for protein coding transcripts, keep track of most extreme start/stop codons, to be used as gene start and gene end
        else if (sep[2].equals("start_codon"))
        {
		    if (strand.equals("+"))
            {
                if (codingStart==0 || codingStart > Integer.parseInt(sep[3]))
                {
                    codingStart = Integer.parseInt(sep[3]);
                }
            }
            else if (codingEnd==0 || codingEnd < Integer.parseInt(sep[4]))
            {
                {
                    codingEnd = Integer.parseInt(sep[4]);
                }
            }
        }
        else if (sep[2].equals("stop_codon"))
        {
		    if (strand.equals("+"))
            {
                if (codingEnd==0 || codingEnd < Integer.parseInt(sep[4]))
                {
                    codingEnd = Integer.parseInt(sep[4]);
                }
            }
            else if (codingStart==0 || codingStart > Integer.parseInt(sep[3]))
            {
                {
                    codingStart = Integer.parseInt(sep[3]);
                }
            }
        }

        //generate splice site intervals using exon start and end positions
        else if (sep[2].equals("exon"))
        {
	        String exonName = sep[8].substring(sep[8].indexOf("exon_id \""));
	        exonName = exonName.substring(9, exonName.indexOf("\";"));
	        int start = Integer.parseInt(sep[3]);
	        int end = Integer.parseInt(sep[4]);

            //update coding start & end
            if (start < geneStart)
            {
                geneStart = start;
            }
            if (end > geneEnd)
            {
                geneEnd = end;
            }

	        //check if duplicate exon
	        duplicateExon = false;
	        for (int i=0; i<donIndex; i++)
            {
	            if (start==donStartPos[i] && end==donEndPos[i])
                {
	                duplicateExon = true;
	            }
	        }
	        for (int i=0; i<accIndex; i++)
            {
	            if (start==accStartPos[i] && end==accEndPos[i])
                {
	                duplicateExon = true;
	            }
	        }
	        if (duplicateExon)
            {
	        	prevLine = "duplicateExon";
	        }
            else
            {
		        //positive strand
		        if (strand.equals("+"))
                {
	                //write donor
	                donStartPos[donIndex]=end-2;
	                donEndPos[donIndex]=end+6;
	                donNames[donIndex]=exonName;
	                donIndex++;

	                //write acceptor
	                if (!(prevLine.equals("gene")||prevLine.equals("transcript")))
                    {
	                    accStartPos[accIndex]=start-20;
	                    accEndPos[accIndex]= start+2;
	                    accNames[accIndex]=exonName;
	                    accIndex++;
	                }
		        }

		        //negative strand
                else if (strand.equals("-"))
                {
	                //write donor
	                donStartPos[donIndex]=start-6;
	                donEndPos[donIndex]=start+2;
	                donNames[donIndex]=exonName;
	                donIndex++;

	                //write acceptor
	                if (!(prevLine.equals("gene")||prevLine.equals("transcript")))
                    {
	                    accStartPos[accIndex]=end-2;
	                    accEndPos[accIndex]=end+20;
	                    accNames[accIndex]=exonName;
	                    accIndex++;
	                }
	            }
		        prevLine = "exon";
	        }
		}
	}

	public static void formatOutputAndReset()
    {
	    //prepare comma delimited strings for accStart, accEnd, donStart and donEnd
	    String accStart = ""; String accEnd = ""; String donStart = ""; String donEnd = ""; String accNameOutput = ""; String donNameOutput = "";
	    boolean multiExon = false;
	    for (int i=0; i<accIndex; i++)
        {
	        multiExon = true;
	        accStart = accStart.concat(Integer.toString(accStartPos[i])).concat(",");
	        accEnd = accEnd.concat(Integer.toString(accEndPos[i])).concat(",");
	        accNameOutput = accNameOutput.concat(accNames[i]).concat(",");
	    }
	    for (int i=0; i<donIndex; i++)
        {
            multiExon = true;
            donStart = donStart.concat(Integer.toString(donStartPos[i])).concat(",");
            donEnd = donEnd.concat(Integer.toString(donEndPos[i])).concat(",");
            donNameOutput = donNameOutput.concat(donNames[i]).concat(",");
        }

        //set gene start and end to coding start and end, if applicable
        if (codingStart!=0)
        {
            geneStart = codingStart;
        }
        if (codingEnd!=0)
        {
            geneEnd = codingEnd;
        }

	    //output splice site intervals for annotated multi-exon genes
	    if (multiExon)
        {
            if (strand.equals("+"))
            {
	    	    outPos[outPosIndex] = chr.concat("\t").concat(Integer.toString(geneStart)).concat("\t").concat(Integer.toString(geneEnd)).concat("\t").concat(geneName).concat("\t").concat("GENE").concat("\t").concat(strand).concat("\t").concat(donStart).concat("\t").concat(donEnd).concat("\t").concat(accStart).concat("\t").concat(accEnd).concat("\t").concat(donNameOutput.concat("\t").concat(accNameOutput));
                outPosIndex++;
            }
            else
            {
	    	    outNeg[outNegIndex] = chr.concat("\t").concat(Integer.toString(geneStart)).concat("\t").concat(Integer.toString(geneEnd)).concat("\t").concat(geneName).concat("\t").concat("GENE").concat("\t").concat(strand).concat("\t").concat(donStart).concat("\t").concat(donEnd).concat("\t").concat(accStart).concat("\t").concat(accEnd).concat("\t").concat(donNameOutput.concat("\t").concat(accNameOutput));
                outNegIndex++;
            }
	    }

	    //reset variables
	    Arrays.fill(accStartPos, 0);
	    Arrays.fill(accEndPos, 0);
	    Arrays.fill(donStartPos, 0);
	    Arrays.fill(donEndPos, 0);
	    Arrays.fill(donNames, "");
	    Arrays.fill(accNames, "");
	    codingStart = 0; codingEnd = 0; donIndex = 0; accIndex = 0;
	}

    public static void writeToFile(String basenameGTF)
    {
        try
        {
            String posFileName = "data/"+basenameGTF+"_SpliceSiteIntervals_pos.txt";
            File file = new File(posFileName);
            FileWriter fw = new FileWriter(file);
            BufferedWriter writer = new BufferedWriter(fw);
            for (int i=0; i<outPosIndex; i++)
            {
                writer.write(outPos[i]+"\n");
            }

            String negFileName = "data/"+basenameGTF+"_SpliceSiteIntervals_neg.txt";
            file = new File(negFileName);
            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            for (int i=0; i<outNegIndex; i++)
            {
                writer.write(outNeg[i]+"\n");
            }

            writer.close();
        }
        catch (IOException e)
        {
            System.out.println(e.getMessage());
        }
    }

}
