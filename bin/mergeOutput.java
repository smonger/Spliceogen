/*
 * Author: Steve Monger 2020
 * 
 * Description: Merges all variant and annotation information to determine motif changes and exon-intron phase.
 * These are integrated to determine logistic regression probabilities.
 * 
 * Input: Stream of gene annotation info, splice site locations, MaxEntScan scores, ESRseq scores, sorted by chr and start position
 * Output: One line per variant with final prediction scores and variant annotation info
 */

import java.io.*;
import java.util.Arrays;
import java.util.jar.Attributes.Name;
import java.lang.*;
import java.text.DecimalFormat;

public class mergeOutput
{
    //store output lines in buffers to minimise I/O
    public static String[] avBuffer = new String[30000];
    public static String[] gainBuffer = new String[30000];
    public static String[] ssBuffer = new String[30000];
    public static int gainBufferIndex = 0;
    public static int avBufferIndex = 0;
    public static int ssBufferIndex = 0;

    // positional distribution scores
    public static double[] donPosBias = new double[550];
    public static double[] accPosBias = new double[550];

    public static void main (String[] args)
    {
        String fileName=args[0];
        String chrAdd=args[1];
        String chrRemove=args[2];
        String strand=args[3];

        //String chrAdd="";
        //String chrRemove="";
        //String fileName="mergeOutTest";
        //String strand="pos";
        writeHeaders(fileName);

        //initialise score tracking variables
        int[] donStart = new int[10000];
        int[] donEnd = new int[10000];
        int[] accStart = new int[10000];
        int[] accEnd = new int[10000];
        String[] donNames = new String[10000];
        String[] accNames = new String[10000];
        double[] scores = new double[24];
        Arrays.fill(scores, -99.0);
        String[] geneID = { "", ".", "", ""};
        String[] prevID = { "", "-99", "", "", "GENE"};
        String s = "";

        //Note array elements:
        //scores[]...mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
        //scores[]...mesDonRef2ndScore(12);mesDonAlt2ndScore(13);mesDonRefPos1(14);mesDonAltPos1(15);mesDonRefPos2(16);mesDonAltPos2(17);
        //scores[]...mesAccRef2ndScore(18);mesAccAlt2ndScore(19);mesAccRefPos1(20);mesAccAltPos1(21);mesAccRefPos2(22);mesAccAltPos2(23);
        //geneID[]...chr(0), name(1), end(2)
        //prevID[]...chr(0), start(1), ref(2), alt(3), type(4)

        //read in positional bias scores
        try
        {
            File file = new File("data/donDist.txt");         
            BufferedReader br = new BufferedReader(new FileReader(file));
            String[] split = new String[2];
            int i=0;
            while ((s = br.readLine()) != null)
            {
                split = s.split("\\s+");
                donPosBias[i] = Double.parseDouble(split[1]);
                i++;
            }	            

            file = new File("data/accDist.txt");         
            br = new BufferedReader(new FileReader(file));
            i=0;
            while ((s = br.readLine()) != null)
            {
                split = s.split("\\s+");
                accPosBias[i] = Double.parseDouble(split[1]);
                i++;
            }
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }	            

        //process input lines
        try
        {
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            while ((s = in.readLine()) != null && s.length() != 0)
            {
                String[] split = s.split("\\s+");
                String chr = split[0];
                int start = Integer.parseInt(split[1]);
                String ref = split[2].toUpperCase();
                String alt = split[3].toUpperCase();
                String type = split[4];

                //check if current line is next variant
                boolean nextVariant = false;
                if (!(chr.equals(prevID[0])&&start==Integer.parseInt(prevID[1])))
                {
                    nextVariant = true;
                }
                else if (!(ref.equals(prevID[2])&&alt.equals(prevID[3])))
                {
                    nextVariant = true;
                }

                //final line flag
                else if (chr.equals("xxx"))
                {
                    nextVariant = true;
                }

                /* Output previous variant if current line is next variant*/

                if(nextVariant && !prevID[4].equals("GENE"))
                {
                    //update overlapping genes and splice sites
                    geneID = updateOverlappingGenes(geneID, prevID);
                    String withinSS = checkForOverlappingSpliceSite(prevID, donNames, accNames, donStart, accStart, donEnd, accEnd, strand);
                    int[] phase = findNearestSpliceSites(prevID, donStart, accStart, donEnd, accEnd, geneID, strand);

                    //format output columns
                    String[] out = new String[35];
                    //scores[]...mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
                    //scores[]...mesDonRef2ndScore(12);mesDonAlt2ndScore(13);mesDonRefPos1(14);mesDonAltPos1(15);mesDonRefPos2(16);mesDonAltPos2(17);
                    //scores[]...mesAccRef2ndScore(18);mesAccAlt2ndScore(19);mesAccRefPos1(20);mesAccAltPos1(21);mesAccRefPos2(22);mesAccAltPos2(23);
                    //out elements
                    //0:chr,1:start,2:end,3:ref,4:alt,5;gene,6:withinSS,7:mesDonRef,8:mesDonAlt,9:mesAccRef,10:mesAccAlt
                    //11:gsDonRef,12:gsDonAlt,13:gsAccRef,14:gsAccAlt,15:ESEmaxRef,16:ESEmaxAlt,17:ESSminRef,18:ESSminAlt
                    //19:donGainP,20:accGainP,21:donLossP,22:accLossP,23:distDon5',24:distDon3',25:distAcc5',26:distAcc3'
                    //27:don1PosRef,28:don1PosAlt,29:don2PosRef,30:don2PosAlt,31:acc1PosRef,32:acc1PosAlt,33:acc2PosRef,34:acc2PosAlt
                    out[0] = prevID[0]; out[1] = Integer.toString(Integer.parseInt(prevID[1]));
                    out[3] = prevID[2]; out[4] = prevID[3]; out[5] = geneID[1]; out[6] = withinSS;

                    //adjust end value for 1bp deletions denoted by "alt=*"
                    int endPos = Integer.parseInt(prevID[1])+prevID[2].length()-1;
                    if (prevID[3].equals("*")) {
                        endPos++;
                    }
                    out[2]= Integer.toString(endPos);

                    //fill scores
                    for (int i=0; i<12; i++)
                    {
                        if (scores[i]==-99 || scores[i]==0)
                        {
                            out[i+7]=".";
                        }
                        else
                        {
                            out[i+7]=Double.toString(scores[i]);
                        }
                    }

                    //update withinSS motif positions
                    scores = updateWithinSSmotifPostions(withinSS, scores);
                    //System.out.println(out[0]+"\t"+out[1]+"\t"+out[3]+"\t"+out[4]+"\t"+strand+"\t"+withinSS+"\t"+Double.toString(scores[20])+"\t"+Double.toString(scores[21])+"\t"+Double.toString(scores[22])+"\t"+Double.toString(scores[23]));

                    //include phase info
                    out[11]= Integer.toString(phase[0]);
                    out[12]= Integer.toString(phase[1]);
                    out[13]= Integer.toString(phase[2]);
                    out[14]= Integer.toString(phase[3]);

                    out[23] = Double.toString(donPhaseScore(phase));
                    out[24] = Double.toString(accPhaseScore(phase));

                    //calculate donor/acceptor creating logistic regression scores
                    String[] lrScores = calculateLogRegScores(scores, out).split("\\s+");
                    out[19]= lrScores[0];
                    out[20]= lrScores[1];
                    out[21]= lrScores[2];
                    out[22]= lrScores[3];

                    /*
                    //temporarily include don/acc score position in output
                    //don
                    for (int i=14; i<18; i++) {
                    if (scores[i]==-99 || scores[i]==0) {
                    out[i+13]=".";
                    } else {
                    out[i+13]=Double.toString(scores[i]);
                    }
                    }
                    //acc
                    for (int i=20; i<24; i++) {
                    if (scores[i]==-99 || scores[i]==0) {
                    out[i+11]=".";
                    } else {
                    out[i+11]=Double.toString(scores[i]);
                    }
                    }
                     */

                    //correct differences in chr nomenclature
                    if (chrAdd.equals("inputAdd=chr"))
                    {
                        out[0]="chr"+out[0];
                    } else if (chrRemove.equals("inputRemove=chr"))
                    {
                        out[0]=out[0].substring(3);
                    }

                    //include delta ESRseq instead of alt scores
                    for (int j=15; j<19; j++)
                    {
                        if(out[j].equals("."))
                        {
                            out[j]="0.00";
                        }
                    }
                    out[15] = String.format("%.02f",Double.parseDouble(out[15])-Double.parseDouble(out[16]));
                    out[16] = String.format("%.02f",Double.parseDouble(out[17])-Double.parseDouble(out[18]));

                    //add line to output buffers and reset scores
                    outputVariantToBuffers(out);
                    Arrays.fill(scores, -99.0);

                    //append to file and empty buffers every ~25000 lines
                    if (avBufferIndex > 25000) {
                        appendToFiles(fileName);
                        resetOutputArrays();
                    }                
                }

                /* Parse Line*/

                //MaxEntScan
                if (type.equals("MESDON"))
                {
                    scores[0]=Double.parseDouble(split[5]);
                    scores[1]=Double.parseDouble(split[6]);
                    scores[18]=Double.parseDouble(split[7]);
                    scores[19]=Double.parseDouble(split[8]);
                    scores[20]=Double.parseDouble(split[9]);
                    scores[21]=Double.parseDouble(split[10]);
                    scores[22]=Double.parseDouble(split[11]);
                    scores[23]=Double.parseDouble(split[12]);
                }
                else if(type.equals("MESACC"))
                {
                    scores[2]=Double.parseDouble(split[5]);
                    scores[3]=Double.parseDouble(split[6]);
                    scores[12]=Double.parseDouble(split[7]);
                    scores[13]=Double.parseDouble(split[8]);
                    scores[14]=Double.parseDouble(split[9]);
                    scores[15]=Double.parseDouble(split[10]);
                    scores[16]=Double.parseDouble(split[11]);
                    scores[17]=Double.parseDouble(split[12]);
                }
                //ESRseq
                else if (type.equals("ESR"))
                {
                    scores[8]=Double.parseDouble(split[5]);
                    scores[9]=Double.parseDouble(split[6]);
                    scores[10]=Double.parseDouble(split[7]);
                    scores[11]=Double.parseDouble(split[8]);
                }
                //GENE
                if (type.equals("GENE"))
                {
                    geneID = updateGeneID(geneID, chr, start, split[3], split[1], split[2]);

                    //if no overlapping genes
                    if (geneID[1].equals(split[3]))
                    {
                        //reset splice site arrays
                        Arrays.fill(donStart, 0);
                        Arrays.fill(donEnd, 0);
                        Arrays.fill(accStart, 0);
                        Arrays.fill(accEnd, 0);
                        Arrays.fill(donNames, null);
                        Arrays.fill(accNames, null);
                    }

                    String[] donStartStr=split[6].split(",");
                    String[] donEndStr=split[7].split(",");
                    String[] accStartStr=split[8].split(",");
                    String[] accEndStr=split[9].split(",");
                    String[] donNamesStr=split[10].split(",");
                    String[] accNamesStr=split[11].split(",");
                    int accIndex = 0;
                    int donIndex = 0;

                    if (split.length>11)
                    {

                        //update splice site pos/name arrays with info from new overlapping gene
                        while (donStart[donIndex]>0)
                        {
                            donIndex++;
                        }
                        while (accStart[accIndex]>0)
                        {
                            accIndex++;
                        }

                        int donMinIndex = Math.min(donStartStr.length, donEndStr.length);
                        int accMinIndex = Math.min(accStartStr.length, accEndStr.length);
                        donMinIndex = Math.min(donMinIndex, donNamesStr.length);
                        accMinIndex = Math.min(accMinIndex, accNamesStr.length);

                        for (int k = 0; k < donMinIndex; k++)
                        {
                            donStart[donIndex+k] = Integer.parseInt(donStartStr[k]);
                            donEnd[donIndex+k] = Integer.parseInt(donEndStr[k]);
                            donNames[donIndex+k] = donNamesStr[k];
                        }
                        for (int k = 0; k < accMinIndex; k++)
                        {
                            accStart[accIndex+k] = Integer.parseInt(accStartStr[k]);
                            accEnd[accIndex+k] = Integer.parseInt(accEndStr[k]);
                            accNames[accIndex+k] = accNamesStr[k];     	            	
                        }
                    }
                }

                //update info for comparison with next line
                prevID[0] = chr;
                prevID[1] = Integer.toString(start);
                prevID[2] = ref;
                prevID[3] = alt;
                prevID[4] = type;
            }

            //final append to file
            if (avBufferIndex > 0)
            {
                appendToFiles(fileName);
            }
            }
            catch (Exception e)
            {
                e.printStackTrace();
            }
        }

        //scans previously overlapping genes, returns updated geneID
        public static String[] updateGeneID(String[] prevGeneID, String chr, int start, String currentGeneName, String currentGeneStart, String currentGeneEnd)
        {
            String[] updatedGeneID = new String[4];

            //first line of input and chromosome changes
            if (prevGeneID[1].isEmpty() || !(prevGeneID[0].equals(chr)))
            {
                updatedGeneID[0] = chr;
                updatedGeneID[1] = currentGeneName;
                updatedGeneID[2] = currentGeneStart;
                updatedGeneID[3] = currentGeneEnd;
                return updatedGeneID;
            }

            String[] geneNameSplit = prevGeneID[1].split(";");
            String[] geneStartSplit = prevGeneID[2].split(";");
            String[] geneEndSplit = prevGeneID[3].split(";");
            updatedGeneID[0] = ""; updatedGeneID[1] = ""; updatedGeneID[2] = ""; updatedGeneID[3] = "";
            int geneIndex = Math.min(geneNameSplit.length, geneEndSplit.length);
            geneIndex = Math.min(geneIndex, geneStartSplit.length);

            for (int i=0; i<geneIndex; i++)
            {
                if (Integer.parseInt(geneEndSplit[i]) >= start)
                {
                    updatedGeneID[1] = updatedGeneID[1].concat(geneNameSplit[i]).concat(";");
                    updatedGeneID[2] = updatedGeneID[2].concat(geneStartSplit[i]).concat(";");
                    updatedGeneID[3] = updatedGeneID[3].concat(geneEndSplit[i]).concat(";");
                }
            }

            //append current gene info
            updatedGeneID[0] = chr;
            updatedGeneID[1] = updatedGeneID[1].concat(currentGeneName);
            updatedGeneID[2] = updatedGeneID[2].concat(currentGeneStart);
            updatedGeneID[3] = updatedGeneID[3].concat(currentGeneEnd);

            return updatedGeneID;
        }

public static String calculateLogRegScores (double[] s, String[] out)
{
    //maxentscan variables
    double imputeVal = -15.0;
    double donGainIntercept = 0.09186; double donGainCoef_c = 0.1266; double donGainCoef_alt = 0.2655;
    double accGainIntercept = -1.0362; double accGainCoef_c = 0.1096; double accGainCoef_alt = 0.2869;
    double donLossIntercept = -0.7309; double donLossCoef_c = -0.7361;
    double accLossIntercept = -0.93825; double accLossCoef_c = -0.8017;
    double donGainWithinIntercept = -0.07592; double donGainWithinCoef_d = 1.3941;
    double accGainWithinIntercept = 0.5292; double accGainWithinCoef_r2 = -0.3682; double accGainWithinCoef_d = 0.5236; double accGainWithinCoef_og = -0.1877;

    //impute missing values
    for (int i=0; i<4; i++)
    {
        if (s[i]==-99.0 | s[i]==0.0)
        {
            s[i] = imputeVal;
        }
    }

    double mesDonChange = s[1] - s[0];
    double mesAccChange = s[3] - s[2];
    double gsDonChange = s[5] - s[4];
    double gsAccChange = s[7] - s[6];
    double pDonGain = -1; double pAccGain = -1; double pDonLoss = -1; double pAccLoss = -1;

    /* Using within/outside splice site logistic regression models, calculate p = 1/(1+e^-(a + b1X1 + b2X2 + ... + bnXn)) */

    // outside splice site
    if (!out[6].contains("ENSE"))
    {
        pDonGain = 1/(1 + Math.exp(-(donGainIntercept + (mesDonChange * donGainCoef_c) + (s[1] * donGainCoef_alt))));
        pAccGain = 1/(1 + Math.exp(-(accGainIntercept + (mesAccChange * accGainCoef_c) + (s[3] * accGainCoef_alt))));

        //phase variable
        double donPhaseIntercept = -8.203503; double donMesLogReg1_coef = 9.941136; double donPhase_coef = 4.040833;
        double accPhaseIntercept = -5.105919; double accMesLogReg1_coef = 6.198231; double accPhase_coef = 7.008072;
        pDonGain = 1/(1 + Math.exp(-(donPhaseIntercept + (pDonGain * donMesLogReg1_coef) + (Double.parseDouble(out[23]) * donGainCoef_alt))));
        pAccGain = 1/(1 + Math.exp(-(accPhaseIntercept + (pAccGain * accMesLogReg1_coef) + (Double.parseDouble(out[24]) * accGainCoef_alt))));
    }

    //within donor site
    if (out[6].contains("donor"))
    {
        pDonLoss = 1/(1 + Math.exp(-(donLossIntercept + (mesDonChange * donLossCoef_c) )));
        double denovoScore = s[13];
        double ogScore = s[3];

        //if highest alt score is at a different site (accounting for indels)
        if ( Math.abs(out[3].length() - out[4].length()) < Math.abs(s[14] - s[15]) )
        {
            denovoScore = s[3];
            ogScore = s[13];
        }
        double denovoChange = denovoScore - s[12];
        pDonGain = 1/(1 + Math.exp(-(donGainWithinIntercept + (denovoScore * donGainWithinCoef_d))));
    }

    //within acceptor site
    if (out[6].contains("acceptor"))
    {
        pAccLoss = 1/(1 + Math.exp(-(accLossIntercept + (mesAccChange * accLossCoef_c) )));
        double denovoScore = s[13];
        double ogScore = s[3];

        //if highest alt score is at a different site (accounting for indels)
        if (Math.abs(out[3].length() - out[4].length()) < Math.abs(s[14] - s[15]) )
        {
            denovoScore = s[3];
            ogScore = s[13];
        }
        double denovoChange = denovoScore - s[12];
        pAccGain = 1/(1 + Math.exp(-(accGainWithinIntercept + ( (s[12] * accGainWithinCoef_r2) + (denovoScore * accGainWithinCoef_d) + (ogScore * accGainWithinCoef_og) ))));
    }

    //format scores
    String pDonGainStr = Double.toString(pDonGain);
    String pAccGainStr = Double.toString(pAccGain);
    String pDonLossStr = Double.toString(pDonLoss);
    String pAccLossStr = Double.toString(pAccLoss);
    try
    {
        pDonGainStr = String.format("%.02f",pDonGain);
        pAccGainStr = String.format("%.02f",pAccGain);
        pDonLossStr = String.format("%.02f",pDonLoss);
        pAccLossStr = String.format("%.02f",pAccLoss);
    }
    catch (Exception e)
    {
        //just use unrounded score
    }
    if (pDonLossStr.equals("-1.00"))
    {
        pDonLossStr = ".";
    }
    if (pAccLossStr.equals("-1.00"))
    {
        pAccLossStr = ".";
    }

    String ret = pDonGainStr + "\t" + pAccGainStr + "\t" + pDonLossStr + "\t" + pAccLossStr;
    return ret;
}

public static void appendToFiles (String fileName) {
    String annovar_file = "temp/"+fileName+"_out_unsorted.txt";
    String gain_file_unsorted = "temp/"+fileName+"_gain_unsorted.txt";
    String ss_file_unsorted = "temp/"+fileName+"_loss_unsorted.txt";

    try {
        //write to _out
        FileWriter fwAV = new FileWriter(annovar_file, true);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        for (int i=0; i<avBufferIndex; i++) {
            avWriter.write(avBuffer[i]+"\n");
        }
        avWriter.close();

        //write to _acceptorCreating
        FileWriter fwGain = new FileWriter(gain_file_unsorted, true);
        BufferedWriter gainWriter = new BufferedWriter(fwGain);
        for (int i=0; i<gainBufferIndex; i++) {
            gainWriter.write(gainBuffer[i]+"\n");
        }
        gainWriter.close();

        //write to _withinSS
        FileWriter fwSS = new FileWriter(ss_file_unsorted, true);
        BufferedWriter ssWriter = new BufferedWriter(fwSS);
        for (int i=0; i<ssBufferIndex; i++)
        {
            ssWriter.write(ssBuffer[i]+"\n");
        }
        ssWriter.close();
    }
    catch (IOException e)
    {
        System.out.println(e.getMessage());
    }
}

public static void writeHeaders (String fileName)
{
    String annovar_file = "output/"+fileName+"_out.txt";
    String gain_file = "temp/"+fileName+"_gain_unsorted.txt";
    String ss_file = "temp/"+fileName+"_loss_unsorted.txt";

    try {
        //write to _out
        FileWriter fwAV = new FileWriter(annovar_file);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        avWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonGainP\taccGainP\tdonLossP\taccLossP\tdistDon5\'\tdistDon3\'\tdistAcc5\'\tdistAcc3\'\tdon1PosRef\tdon1PosAlt\tdon2PosRef\tdon2PosAlt\tacc1PosRef\tacc1PosAlt\tacc2PosRef\tacc2PosAlt"+"\n");
        avWriter.close();
    }
    catch (IOException e)
    {
        System.out.println(e.getMessage());
    }
} 

public static void resetOutputArrays()
{
    avBuffer = new String[30000];
    gainBuffer = new String[30000];
    ssBuffer = new String[30000];
    avBufferIndex = 0; gainBufferIndex = 0; ssBufferIndex = 0;
    System.gc();
}

//if previous genes overlapped, check individually whether they overlap this variant and update
public static String[] updateOverlappingGenes(String[] geneID, String[] prevID)
{
    if (geneID[1].contains(";"))
    {
        String[] geneID_original = geneID;
        String[] geneNameSplit = geneID[1].split(";");
        String[] geneStartSplit = geneID[2].split(";");
        String[] geneEndSplit = geneID[3].split(";");
        geneID[1]="";
        geneID[2]="";
        geneID[3]="";
        int geneIndex = Math.min(geneEndSplit.length, geneNameSplit.length);
        geneIndex = Math.min(geneIndex, geneStartSplit.length);

        for (int i=0; i<geneIndex; i++)
        {
            if (Integer.parseInt(geneEndSplit[i]) >= Integer.parseInt(prevID[1]))
            {
                if (!geneID[1].equals(""))
                {
                    geneID[1] = geneID[1].concat(";");
                    geneID[2] = geneID[2].concat(";");
                    geneID[3] = geneID[3].concat(";");
                }
                geneID[1] = geneID[1].concat(geneNameSplit[i]);
                geneID[2] = geneID[2].concat(geneStartSplit[i]);
                geneID[3] = geneID[3].concat(geneEndSplit[i]);
            }
        }
        if (geneID[1].equals(""))
        {
            return geneID_original;
        }
    }
    return geneID;
}

public static String checkForOverlappingSpliceSite( String[] prevID, String[] donNames, String[] accNames, int[] donStart, int[] accStart, int[] donEnd, int[] accEnd, String strand)
{
    String withinSS="";
    int startPos = Integer.parseInt(prevID[1]);
    int endPos = Integer.parseInt(prevID[1])+prevID[2].length()-1;
    if (prevID[3].equals("*")) {
        endPos++;
    }
    for (int i=0; i<donNames.length; i++)
    {
        if (donNames[i]==null)
        {
            break;
        }
        if (!withinSS.contains(donNames[i]))
        {
            if (endPos>=donStart[i]&&startPos<=donEnd[i])
            {
                if (!withinSS.equals(""))
                {
                    withinSS=withinSS.concat(",");
                }
                withinSS = withinSS.concat(donNames[i]).concat("_donor");

                //identify postion in motif
                String mPos = ""; int motifPos = 0;
                if (strand.equals("pos"))
                {
                    motifPos = startPos - donStart[i] - 3;
                }
                else
                {
                    motifPos = donEnd[i] - startPos - 3;
                }
                if (motifPos >= 0)
                {
                    motifPos++;
                    mPos="I"+Integer.toString(motifPos);
                }
                else
                {
                    motifPos = motifPos * -1;
                    mPos="E"+Integer.toString(motifPos);
                }
                withinSS = withinSS.concat("("+mPos+")");
            }
        }
    }

    for (int i=0; i<accNames.length; i++)
    {
        if (accNames[i]==null)
        {
            break;
        }
        if (!withinSS.contains(accNames[i]))
        {
            if (endPos>=accStart[i]&&startPos<=accEnd[i])
            {
                if (!withinSS.equals("")) {
                    withinSS=withinSS.concat(",");
                }
                withinSS = withinSS.concat(accNames[i]).concat("_acceptor");

                //identify postion in motif
                String mPos = ""; int motifPos = 0;
                if (strand.equals("pos"))
                {
                    motifPos = accEnd[i] - startPos - 3;
                }
                else
                {
                    motifPos = startPos - accStart[i] - 3;
                }
                if (motifPos >= 0)
                {
                    motifPos++;
                    mPos="I"+Integer.toString(motifPos);
                }
                else
                {
                    motifPos = motifPos * -1;
                    mPos="E"+Integer.toString(motifPos);
                }
                withinSS = withinSS.concat("("+mPos+")");
            }
        }
    }
    if (withinSS.equals(""))
    {
        withinSS=".";
    }
    return withinSS;
}

// outputs distance to 4 neighbouring splice sites (5' donor, 3' donor, 5' acceptor, 3' acceptor)
// if within first/last exon, record distance to transcript start/end instead. Negative
public static int[] findNearestSpliceSites( String[] prevID, int[] donStart, int[] accStart, int[] donEnd, int[] accEnd, String[] geneID, String strand)
{
    int startPos = Integer.parseInt(prevID[1]);
    int[] phase = new int[4];
    int i = 0;
    int geneStart = getMinGeneStart(geneID);
    int geneEnd = getMaxGeneEnd(geneID);

    //positive strand
    if (strand.equals("pos"))
    {
        phase[0] = -9999999; phase[2] = -9999999;
        phase[1] = Integer.MAX_VALUE; phase[3] = Integer.MAX_VALUE;

        while (donStart[i] > 0)
        {
            //closest upstream donor
            if (donStart[i]+2 <= startPos && phase[0] < donStart[i]+2)
            {
                phase[0] = donStart[i]+2;
            }

            //closest downstream donor
            else if (donStart[i]+2 > startPos && phase[1] > donStart[i]+2)
            {
                phase[1] = donStart[i]+2;	                	
            }
            i++;
        }

        i = 0;
        while (accStart[i] > 0)
        {
            //closest upstream acceptor
            if (accStart[i]+20 <= startPos && phase[2] < accStart[i]+20)
            {
                phase[2] = accStart[i]+20;
            }

            //closest downstream acceptor
            else if (accStart[i]+20 > startPos && phase[3] > accStart[i]+20)
            {
                phase[3] = accStart[i]+20;
            }
            i++;
        }

        //if within first exon/intron (no 5' acceptor)
        if (phase[2] == -9999999)
        {
            //check for conflicting overlapping gene end positions
            phase[2] = geneStart - startPos;

            //if within first exon (no 5' donor)
            if (phase[0] == -9999999)
            {
                phase[0] = geneStart - startPos;
            }
			else
			{
				phase[0] = startPos - phase[0];
			}
        }
		else
		{
            phase[0] = startPos - phase[0];
			phase[2] = startPos - phase[2];
		}

        //if within last exon/intron (no 3' donor)
        if (phase[1] == Integer.MAX_VALUE)
        {
            //check for conflicting overlapping gene end positions
            phase[1] = startPos - geneEnd;

            //if within last exon (no 3' acceptor)
            if (phase[3] == Integer.MAX_VALUE)
            {
                phase[3] = startPos - geneEnd;
            }
			else
			{
				phase[3] = phase[3] - startPos;
			}
        }
		else
		{
			phase[1] = phase[1] - startPos;
            phase[3] = phase[3] - startPos;
		}
    }

    //negative strand
    else if (strand.equals("neg"))
    {
        phase[0] = Integer.MAX_VALUE; phase[2] = Integer.MAX_VALUE;
        phase[1] = -9999999; phase[3] = -9999999;

        while (donStart[i] > 0)
        {
            //closest upstream donor
            if (donStart[i]+6 > startPos && phase[0] > donStart[i]+6)
            {
                phase[0] = donStart[i]+6;
            }
            //closest downstream donor
            else if (donStart[i]+6 <= startPos && phase[1] < donStart[i]+6)
            {
                phase[1] = donStart[i]+6;	                	
            }
            i++;
        }

        i = 0;
        while (accStart[i] > 0)
        {
            //closest upstream acceptor
            if (accStart[i]+2 >= startPos && phase[2] > accStart[i]+2)
            {
                phase[2] = accStart[i]+2;
            }

            //closest downstream acceptor
            else if (accStart[i]+2 < startPos && phase[3] < accStart[i]+2)
            {
                phase[3] = accStart[i]+2;
            }
            i++;
        }

        //if within first exon/intron (no 5' acceptor) 
        if (phase[2] == Integer.MAX_VALUE)
        {		        	
            phase[2] = startPos - geneEnd;

            //if within first exon (no 5' donor)
            if (phase[0] == Integer.MAX_VALUE)
			{
                phase[0] = startPos - geneEnd;
            }
			else
			{
				phase[0] = phase[0] - startPos ;
			}
        }
		else
		{
			phase[0] = phase[0] - startPos ;
			phase[2] = phase[2] - startPos;
		}

        //if within last exon/intron (no 3' donor) 
        if (phase[1] == -9999999)
        {
            phase[1] = geneStart - startPos;

            //if within last exon (no 3' acceptor)
            if (phase[3] == -9999999)
            {
                phase[3] = geneStart - startPos;
            }
			else
			{
				phase[3] = startPos - phase[3];
			}
        }      
		else
		{
			phase[1] = startPos - phase[1];
			phase[3] = startPos - phase[3];
		}
    }   
    return phase;
}	    

public static void outputVariantToBuffers(String[] out) {	
    //write full line for annovar
    String avLine = "";
    for (int i=0; i<24; i++)
    {
        avLine = avLine+out[i]+"\t";
    }
    avLine = avLine+out[24];
    avBuffer[avBufferIndex] = avLine;
    avBufferIndex++;
    /*	        //write withinSS
                if (out[6].contains("ENSE")) {
                String ssLine = out[0]+"\t"+out[1]+"\t"+out[3]+"\t"+out[4]+"\t"+out[5]+"\t"+out[6]+"\t"+out[21]+"\t"+out[22];
                Double lossMax = 0.00 ;
                try {
                if (!out[21].equals(".")) {
                lossMax = Double.parseDouble(out[21]);
                }
                if (!out[22].equals(".")) {
                if (lossMax < Double.parseDouble(out[22])) {
                lossMax = Double.parseDouble(out[22]);
                }
                }
                } catch (Exception e) {}
                ssLine = ssLine+"\t"+lossMax;
                ssBuffer[ssBufferIndex] = ssLine;
                ssBufferIndex++;
                } else {
    //write gain
    String gainLine = out[0]+"\t"+out[1]+"\t"+out[3]+"\t"+out[4]+"\t"+out[5]+"\t"+out[19]+"\t"+out[20];
    Double gainMax = 0.00 ;
    try {
    if (!out[19].equals(".")) {
    gainMax = Double.parseDouble(out[19]);
    }
    if (!out[20].equals(".")) {
    if (gainMax < Double.parseDouble(out[20])) {
    gainMax = Double.parseDouble(out[20]);
    }
    }
    }
    catch (Exception e) {}
    if (gainMax>=0.7) {
    gainLine = gainLine+"\t"+Double.toString(gainMax);
    gainBuffer[gainBufferIndex] = gainLine;
    gainBufferIndex++;
    }
    }              
     */	    }

    // returns "maximum gene end" int value when overlapping genes may have distinct end positions
    // 
public static int getMaxGeneEnd(String[] geneID)
{	    	
    int pos = 0;
    if (geneID[3].contains(";"))
    {
        String[] split = geneID[3].split(";");
        pos = Integer.parseInt(split[0]);
        for (int j=1; j<split.length; j++)
        {
            if (Integer.parseInt(split[j]) < pos)
            {
                pos = Integer.parseInt(split[j]);
            }
        }
    }
    else if (!geneID[3].equals(""))
    {
        pos = Integer.parseInt(geneID[3]);
    }	    	
    return pos;
}

// returns "minimum gene start" int value when overlapping genes may have distinct start positions
// 	    
public static int getMinGeneStart(String[] geneID)
{	    	
    int pos = 0;
    if (geneID[2].contains(";"))
    {
        String[] split = geneID[2].split(";");
        pos = Integer.parseInt(split[0]);
        for (int j=1; j<split.length; j++)
        {
            if (Integer.parseInt(split[j]) < pos)
            {
                pos = Integer.parseInt(split[j]);
            }
        }
    }
    else if (!geneID[3].equals(""))
    {
        pos = Integer.parseInt(geneID[2]);
    }	    	
    return pos;
}

public static double[] updateWithinSSmotifPostions( String withinSS, double[] scores)
{
    //withinSSgain get position of new site
    if (withinSS.contains("donor"))
    {
        String varPosStr = withinSS.substring(withinSS.indexOf("d") + 1);
        varPosStr = varPosStr.substring(varPosStr.indexOf("(") + 1);
        varPosStr = varPosStr.substring(0, varPosStr.indexOf(")"));
        int varPos = Integer.parseInt(String.valueOf(varPosStr.charAt(1)));

        if (varPosStr.charAt(0) == 'E')
        {
            for (int j = 0; j< 4; j++)
            {
                if (scores[j+20] != -99)
                {
                    scores[j+20] = scores[j+20] + varPos - 9;
                    if (scores[j+20] <= 0)
                    {
                        scores[j+20] = -1 - scores[j+20];
                    }
                }
            }
        } else if (varPosStr.charAt(0) == 'I') {
            for (int j = 0; j< 4; j++) {
                if (scores[j+20] != -99) {
                    scores[j+20] = (varPos + 9 - scores[j+20]) * -1;
                    if (scores[j+20] >= 0) {
                        scores[j+20] = scores[j+20] + 1;
                    }
                }
            }
        }
    }
    if (withinSS.contains("acceptor"))
    {
        String varPosStr = withinSS.substring(withinSS.indexOf("a") + 1);
        varPosStr = varPosStr.substring(varPosStr.indexOf("(") + 1);
        varPosStr = varPosStr.substring(0, varPosStr.indexOf(")"));
        int varPos = Integer.parseInt(String.valueOf(varPosStr.charAt(1)));

        if (varPosStr.charAt(0) == 'E')
        {
            for (int j = 0; j< 4; j++)
            {
                if (scores[j+20] != -99)
                {
                    scores[j+20] = scores[j+20] + varPos - 20;
                    if (scores[j+20] <= 0)
                    {
                        scores[j+20] = -1 - scores[j+20];
                    }
                }
            }
        }
        else if (varPosStr.charAt(0) == 'I')
        {
            for (int j = 0; j< 4; j++)
            {
                if (scores[j+20] != -99)
                {
                    scores[j+20] = (varPos + 20 - scores[j+20]) * -1;
                    if (scores[j+20] >= 0)
                    {
                        scores[j+20] = scores[j+20] + 1;
                    }
                }
            }
        }
    }
    return scores;
}

// input: distances to 4 neighbouring splice sites
// returns a score representing likelihood of donor gain based on phase
//
public static double donPhaseScore(int[] phase)
{
    double val = 0;

    //within final exon or final intron trough
    if (phase[1] < 0 && phase[3] < 80)
    {
        val = 0.04;
    }

    //within any other trough
    else if (phase[3] > 0 && phase[3] < 80)
    {
        val = donPosBias[200-phase[3]];
    }

    //within first exon/intron peak
    else if (phase[2] < 0 && phase[2] > -350)
    {
        val = donPosBias[200-phase[2]];
    }

    //within internal peak
    else if (phase[2] > 0 && phase[2] < 350)
    {
        val = donPosBias[200+phase[2]];
    }

    //within secondary peak
    else if (phase[3] > 0 && phase[3] < 200)
    {
        val = donPosBias[200-phase[3]];
    }

    //intronic baseline
    else
    {
        val = 0.13;
    }

    return val;
}

// input: distances to 4 neighbouring splice sites
// returns a score representing likelihood of acceptor gain based on phase
//
public static double accPhaseScore(int[] phase)
{
    double val = 0;

    //within first exon or first intron trough
    if (phase[2] < 0 && phase[0] < 60)
    {
        val = 0.02;
    }

    //within any other trough
    else if (phase[0] < 80)
    {
        val = accPosBias[350+phase[0]];
    }

    //within last exon/intron peak
    else if (phase[1] < 0 && phase[1] > -350)
    {
        val = accPosBias[350+phase[1]];
    }

    //within internal peak
    else if (phase[1] < 350 && phase[1] > 0)
    {
        val = accPosBias[350-phase[1]];
    }

    //within secondary peak
    else if (phase[0] >= 80 && phase[0] < 200)
    {
        val = accPosBias[350+phase[0]];
    }

    //intronic baseline
    else
    {
        val = 0.09;
    }

    return val;
}

}
