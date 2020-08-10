#!/bin/bash

function printVersion {
    printf "Spliceogen 2.0 1-October-2019\n"
}

function printHelp {
    cat <<-END
Usage:
------
3 required args:
1)  -input path/to/VCF/input/file(s).VCF
        Note: multiple input files are accepted "eg. -input *.vcf"
        Note: deprecated v1.0 tags "inputVCF" and "inputBED" are still accepted
2)  -gtf path/to/annotation.gtf
3)  -fasta path/to/genome.fa
optional arg:
4)  -branchpointer hgXX
        OR
    -branchpointerIndels hgXX
        Note: user must specify hg19 or hg38
END
}

#set default parameters
POSITIONAL=()
INPUTFILES=""
USEBP=""
USEBPINDELS=""
ANNOTATION=""
FASTAPATH=""
GENOMEBUILD=""
STRANDS="pos neg"
#parse command line args
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    printHelp
    exit 0
    ;;
    -v|--version)
    printVersion
    exit 0
    ;;
    -input)
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -inputVCF)
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -inputBED)
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -fasta)
    FASTAPATH="$2"
    shift
    shift ;;
    -gtf)
    ANNOTATION="$2"
    shift
    shift ;;
    -branchpointer)
    GENOMEBUILD="$2"
    USEBP="TRUE"
    USEBPINDELS="FALSE"
    shift
    shift ;;
    -branchpointerIndels)
    GENOMEBUILD="$2"
    USEBP="TRUE"
    USEBPINDELS="TRUE"
    shift
    shift ;;
    -strand)
    STRANDS="$2"
    shift
    shift ;;
    *)
    POSITIONAL+=("$1")
    shift ;;
esac
done
set - "${POSITIONAL[@]}"

#check input files exist
if [ ! -f $FASTAPATH ]; then
    echo -e "Fasta file not found: use -fasta ./path/to/hgXX.fa\nExiting..."
    exit 1
elif [ ! -f "$ANNOTATION" ]; then
    echo "GTF annotation file not found: use -gtf path/to/gencodeXX.gtf\nExiting..."
    exit 1
fi

#correct mismatches in "chr" nomenclature between gtf and fasta
gtfChr=$(zcat -f "$ANNOTATION" | grep -v '^GL000' | tail -1 | awk '{print $1}' | grep chr)
fastaChr=$(head -1 "$FASTAPATH" | awk '{print $1}' | grep chr)
    gtfChrAdd=""
    gtfChrRemove="UnmatchedString"
    if [ "$fastaChr" != "" ]; then
        if [ "$gtfChr" == "" ]; then
            gtfChrAdd="chr"
        fi        
    elif [ "$fastaChr" == "" ]; then
        if [ "$gtfChr" != "" ]; then
            gtfChrRemove="chr"
        fi        
    fi

#prepare splice site intervals from annotation
gtfBasename=$(basename $ANNOTATION)
if [ ! -f data/"$gtfBasename"_SpliceSiteIntervals_pos.txt ] || [[ "$ANNOTATION" -nt data/"$gtfBasename"_SpliceSiteIntervals_pos.txt ]] ; then
    echo "Preparing splice site annotation..."
    zcat -f "$ANNOTATION" | grep -v '^GL000' | sed "s/$gtfChrRemove//" | awk -v var="$gtfChrAdd" '{print var$0}' | java -cp bin getSpliceSiteIntervalsFromGTF "$gtfBasename"
    #zcat -f "$ANNOTATION" | grep -v '^GL000' | sed "s/$gtfChrRemove//" | awk -v var="$gtfChrAdd" '{print var$0}' | java -cp bin getSpliceSiteIntervalsFromGTF | grep -E "[[:space:]]\-[[:space:]]" > data/"$gtfBasename"_SpliceSiteIntervals_neg.txt
fi

#for each input VCF/TSV file
for FILE in $INPUTFILES; do
    fileID=$(echo "$FILE" | xargs -n 1 basename)

    #check current file exists 
    if [ ! -f "$FILE" ]; then
        echo "Error: variant input file not found: $FILE \n Exiting..."
        exit 1
    fi
    echo "Input file: $fileID"

    #write headers
    echo -e "#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonGainP\taccGainP\tdonLossP\taccLossP" > output/"$fileID"_out.txt
    echo -e "#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonGainP\taccGainP\tdonLossP\taccLossP\tdistDon5\'\tdistDon3\'\tdistAcc5\'\tdistAcc3\'\tdon1PosRef\tdon1PosAlt\tdon2PosRef\tdon2PosAlt\tacc1PosRef\tacc1PosAlt\tacc2PosRef\tacc2PosAlt" > output/"$fileID"_out.txt

    #remove temp files from any previous run
    rm temp/"$fileID"* output/"$fileID"mesOmmitted.txt output/"$fileID"refMismatch.txt 2> /dev/null

    #determine input file type
    FILETYPE=""
    nFields=$(tail -1 < <(zcat -f "$FILE"))
    nFields=$(echo "$nFields" | wc -w)
    #vcfHeader=$(zcat -f $FILE | head -1 | grep VCF)
    if [ "$nFields" -eq 4 ]; then
        FILETYPE="TSV"
    else
    #elif [ ! -z "$vcfHeader" ]; then
        FILETYPE="VCF"
    fi

    #correct mismatches in "chr" nomenclature between variant input and provided gtf/fasta
    #inputChr=$(zcat -f "$FILE" | tail -1 | awk '{print $1}' | grep chr)
    inputChr=$(tail -1 < <(zcat -f "$FILE"))
    inputChr=$(echo "$inputChr" | awk '{print $1}' | grep chr)
    inputChrAdd=""
    inputChrRemove="UnmatchedString"
    if [ "$fastaChr" != "" ]; then
        if [ "$inputChr" == "" ]; then
            inputChrAdd="chr"
        fi        
    elif [ "$fastaChr" == "" ]; then
        if [ "$inputChr" != "" ]; then
            inputChrRemove="chr"
        fi        
    fi

    #sort body of input file
    zcat -f "$FILE" | grep "^#" > temp/"$fileID"_sorted
    if [ "$FILETYPE" == "TSV" ]; then
        zcat -f "$FILE" | awk NF | grep -v "^#" | sort -k1,1 -k2,2n | sed "s/$inputChrRemove//" | awk -v OFS="\\t" -v var=$inputChrAdd '{print var$1, $2, $2, "x", "1", ".", $3, $4}' >> temp/"$fileID"_sorted
    else 
        zcat -f "$FILE" | grep -v "^#" | sort -k1,1 -k2,2n | sed "s/$inputChrRemove//" | awk -v OFS="\\t" -v var=$inputChrAdd '{print var$0}' >> temp/"$fileID"_sorted
    fi

    #check bedtools is installed
    bedtoolsLocation=$(which bedtools);
    if [ "$bedtoolsLocation" == "" ]; then
        printf -- 'Warning: Bedtools does not appear to be installed.\n';
        printf -- 'Get it here: https://bedtools.readthedocs.io/en/latest/content/installation.html\n';
    fi;

    #note: the bedtools getfasta "-name" behaviour is not backwards compatible. Between v2.26.0 and v2.27.0 the previous function  of "name" was given to "name+". So need to alter command for different versions
    recentBedtools=false
    bedtoolsVersion=$(bedtools -version)
    versionSort=$(echo -e "bedtools v2.26.0\nbedtools v2.27.0\n$bedtoolsVersion" | sort -V | head -2 | tail -1)
    if [ "$versionSort" == "bedtools v2.27.0" ]; then
        recentBedtools=true
    fi

    for strand in $STRANDS; do
        strandSym="+"
        if [ "$strand" == "neg" ]; then
            strandSym="-"
        fi
        echo "Running "$strand" strand..."

        #bedtools intersect to exclude intergenic variants
        zcat -f "$ANNOTATION"| awk -v OFS="\\t" -v var="$strandSym" '$7 == var && $3 == "gene" {print}' | sort -k1,1 -k4,4n | grep -v '^GL000' | awk -v var="$gtfChrAdd" -v OFS="\\t" '{print var$0}' | sed "s/$gtfChrRemove//" | bedtools intersect -a temp/"$fileID"_sorted -b stdin -wa -wb -sorted  > temp/"$fileID"unstrandedInput.txt 
        if [ $? -ne 0 ]; then
            echo "Bedtools intersect returned non-zero exit status. Intersection failed between provided variant VCF/BED file and provided GTF. Exiting..."
            exit 1
        fi

        #generate flanking intervals.bed for bedtools getfasta
        if [ "$FILETYPE" = "VCF" ]; then
            awk -v var="$strandSym" -v OFS="\\t" '{print ".", $1, $2, var, $4, $5}' temp/"$fileID"unstrandedInput.txt | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        elif [ "$FILETYPE" = "TSV" ]; then
            awk -v var="$strandSym" -v OFS="\\t" '{print ".", $1, $2, var, $7, $8}' temp/"$fileID"unstrandedInput.txt | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        fi

        echo "Retrieving flanking FASTA sequence..."
        if [ "$recentBedtools" == true ]; then
            bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name+ -s > temp/"$fileID"seqToScan.FASTA
        else
            bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name -s > temp/"$fileID"seqToScan.FASTA
        fi
        if [ ! -s temp/"$fileID"seqToScan.FASTA ]; then
            echo "Warning: no variants were returned following bedtools getfasta command on $strand strand"
        fi

        #seqScan: generates input strings for maxentscan and ESRseq scores
        echo "Scanning for motifs..."
        java -cp bin seqScan temp/"$fileID"seqToScan.FASTA -useESR $fileID >&2

        #run maxEntScan and confirm non-zero exit, since invalid inputs cause it to exit early
        if [ -s temp/"$fileID"mesDonorInput.txt ] || [ -s temp/"$fileID"mesAcceptorInput.txt ] ; then
            echo "Running MaxEntScan..."
            perl score5.pl temp/"$fileID"mesDonorInput.txt | sed 's/;-;(-)/;/' | sed 's/;+;/;/' | sort -t '>' -k2 | java -cp bin processScoresMES > temp/"$fileID"mesDonorScores.txt
            retVal=( ${PIPESTATUS[0]} )
            if [ $retVal -ne 0 ]; then
                echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
            exit $retVal
            fi
            perl score3.pl temp/"$fileID"mesAcceptorInput.txt | sed 's/;-;(-)/;/' | sed 's/;+;/;/' | sort -t '>' -k2 | java -cp bin processScoresMES > temp/"$fileID"mesAcceptorScores.txt
            retVal=( ${PIPESTATUS[0]} )
            if [ $retVal -ne 0 ]; then
                echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
            exit $retVal
            fi
        fi

        #determine files to merge
        scoresToMerge=""
        if [ -s temp/"$fileID"mesDonorScores.txt ] ; then
            scoresToMerge="$scoresToMerge temp/"$fileID"mesDonorScores.txt"
        fi
        if [ -s temp/"$fileID"mesAcceptorScores.txt ] ; then
            scoresToMerge="$scoresToMerge temp/"$fileID"mesAcceptorScores.txt"
        fi
        if [ -s temp/"$fileID"ESRoutput.txt ] ; then
            scoresToMerge="$scoresToMerge temp/"$fileID"ESRoutput.txt"
        fi

        #edit splice site intervals file in event of changed input/fasta "chr" nomenclature
        intervalsFileChr=$(head -1 data/"$gtfBasename"_SpliceSiteIntervals_"$strand".txt | awk '{print $1}' | grep chr)
        if [ "$fastaChr" != "" ]; then
            if [ "$intervalsFileChr" == "" ]; then
                sed -i 's/^/chr/' data/"$gtfBasename"_SpliceSiteIntervals_"$strand".txt
            fi
        elif [ "$fastaChr" == "" ]; then
            if [ "$gtfChr" != "" ]; then
                sed -i "s/^chr//" data/"$gtfBasename"_SpliceSiteIntervals_"$strand".txt
            fi
        fi
        checkScoresExist=$(echo "$scoresToMerge" | grep "temp")
        if [ -z "$checkScoresExist" ]; then
            echo "No MaxEntScan/GeneSplicer/ESRseq scores to process"
        else
            echo "Processing scores..."
            cat $(echo "$scoresToMerge") data/"$gtfBasename"_SpliceSiteIntervals_"$strand".txt sources/terminatingMergeLine.txt | sort -k1,1 -V -k 2,2n -k 3 -k 4 -s | java -cp bin mergeOutput "$fileID" inputAdd="$inputChrRemove" inputRemove="$inputChrAdd" "$strand"

        fi
	    rm temp/"$fileID"mes* temp/"$fileID"ESRoutput.txt
    done

    #report FASTA issues- ref allele mismatches and invalid characters
    if [ -s output/"$fileID"refMismatch.txt ]; then
        refMismatchCount=$(sort output/"$fileID"refMismatch.txt | sed 's/-;(-)/+;/g' | uniq | wc -l | awk '{print $1}')
        echo "Note: $refMismatchCount variant(s) were excluded because the provided Reference allele does not match the nucleotides in the provided FASTA. IDs of excluded variants are listed here: Spliceogen/output/""$fileID""refMismatch.txt"  >&2
    fi
    if [ -s output/"$fileID"mesOmmitted.txt ]; then
        mesOmmittedCount=$(sort output/"$fileID"mesOmmitted.txt | sed 's/-;(-)/+;/g' | uniq | wc -l | awk '{print $1}')
        echo "Note: $mesOmmittedCount variant(s) were excluded because their flanking FASTA sequence contains invalid characters (most commonly \"n\"). IDs of ommitted variants are listed here: Spliceogen/output/""$fileID""mesOmmitted.txt"
    fi

    #sort predictions
    if [ -s temp/"$fileID"_out_unsorted.txt ]; then
        echo "sorting predictions..."
        #set if need ESR scores
        if [ "useESR" == "true" ]; then
            echo -e "#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tdistDon5'\tdistDon3'\tdistAcc5'\tdistAcc3'\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonGainP\taccGainP\tdonLossP\taccLossP\tdonPhaseScore\taccPhaseScore" > output/"$fileID"_out.txt
            sort -k1,1 -V -k 2,2n -k 3 -k 4 -s output/"$fileID"_unsortedBothStrands.txt >> output/"$fileID"_out.txt
        else
            echo -e "#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tdistDon5'\tdistDon3'\tdistAcc5'\tdistAcc3'\tdESE\tdESS\tdonGainP\taccGainP\tdonLossP\taccLossP\tdonPhaseScore\taccPhaseScore" > output/"$fileID"_out.txt
            sort -k1,1 -V -k 2,2n -k 3 -k 4 -s temp/"$fileID"_out_unsorted.txt | cut --complement -f3,8,9,10,11,18,19 >> output/"$fileID"_out.txt
        fi
    fi
    #clean up temp files
    rm temp/"$fileID"* 2> /dev/null
    #cat output/"$fileID"_out.txt | cut --complement -f11-19 > output/"$fileID"_out.txt_temp
    #rm output/"$fileID"_out.txt

    #echo -e "#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tdonGainP\taccGainP\tdonLossP\taccLossP\tdistDon5\'\tdistDon3\'\tdistAcc5\'\tdistAcc3\'\tdon1PosRef\tdon1PosAlt\tdon2PosRef\tdon2PosAlt\tacc1PosRef\tacc1PosAlt\tacc2PosRef\tacc2PosAlt" > output/"$fileID"_out.txt
    #grep -v "^#" output/"$fileID"_out.txt_temp | awk -v OFS="\\t" '{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14}' >> output/"$fileID"_out.txt
    #rm output/"$fileID"_unsortedBothStrands.txt output/*ssGain* output/*withinSS*


#start comment to skip duplicate operations
#    #merge predictions from both strands
#    while read -r chr start end ref alt gene within donGain accGain donLoss accLoss ; do
#	#if not first line
#        if [ "$prevLine" != "" ]; then
#            newDonGain=$(echo "$donGain")
#            newAccGain=$(echo "$accGain")
#            prevStart=$(echo "$prevLine" | awk '{print $2}')
#	    duplicate="false"
#	    #check for duplicate
#            if [ "$prevStart" == "$start" ]; then
#            	prevRef=$(echo "$prevLine" | awk '{print $4}')
#                if [ "$prevRef" == "$ref" ]; then
#            	    prevAlt=$(echo "$prevLine" | awk '{print $5}')
#                    if [ "$prevAlt" == "$alt" ]; then
#		        duplicate="true"
#		        fi
#		    fi
#	    fi
#
#	    #if duplicate line
#	    if [ "$duplicate" == "true" ]; then
#                #overlapping, so find highest scores
#                prevGene=$(echo "$prevLine" | awk '{print $6}')
#		allGenes=$(echo "$prevGene,$gene")
#                #donGain
#                oldDonGain=$(echo "$prevLine" | awk '{print $8}')
#		if [ 1 -eq "$(echo "${donGain} < ${oldDonGain}" | bc)" ]; then
#                    newDonGain=$(echo "$oldDonGain")
#		fi
#                #accGain
#                oldAccGain=$(echo "$prevLine" | awk '{print $9}')
#		if [ 1 -eq "$(echo "${accGain} < ${oldAccGain}" | bc)" ]; then
#                    newAccGain=$(echo "$oldAccGain")
#		fi
#		#withinSS
#                oldWithin=$(echo "$prevLine" | awk '{print $7}')
#		newWithin=$(echo "$oldWithin,$within")
#		newDonLoss=$(echo "$donLoss")
#		newAccLoss=$(echo "$accLoss")
#                oldDonLoss=$(echo "$prevLine" | awk '{print $10}')
#                oldAccLoss=$(echo "$prevLine" | awk '{print $11}')
#		#donLoss
#		if [ "$oldDonLoss" != "." ]; then
#		    if [ "$newDonLoss" != "." ]; then
#		        #scores for both strands so need to compare
#			if [ 1 -eq "$(echo "${donLoss} < ${oldDonLoss}" | bc)" ]; then
#			    newDonLoss=$(echo "$oldDonLoss")
#			fi
#		    else
#			newDonLoss=$(echo "$oldDonLoss")
#		    fi
#		fi
#		#accLoss
#		if [ "$oldAccLoss" != "." ]; then
#		    if [ "$newAccLoss" != "." ]; then
#		        #scores for both strands so need to compare
#			if [ 1 -eq "$(echo "${accLoss} < ${oldAccLoss}" | bc)" ]; then
#			    newAccLoss=$(echo "$oldAccLoss")
#			fi
#		    else
#			newAccLoss=$(echo "$oldAccLoss")
#		    fi
#		fi
#		#set updated values
#	    	prevLine=$(echo -e "$chr\t$start\t$end\t$ref\t$alt\t$allGenes\t$newWithin\t$newDonGain\t$newAccGain\t$newDonLoss\t$newAccLoss")
#	    #not duplicate
#	    else
#    		echo -e "$prevLine" >> output/"$fileID"_out.txt_merged
#	    	prevLine=$(echo -e "$chr\t$start\t$end\t$ref\t$alt\t$gene\t$within\t$donGain\t$accGain\t$donLoss\t$accLoss")
#	    fi
#
#	#first line only
#	else
#	    rm "$file"_fixed
#	    prevLine=$(echo -e "$chr\t$start\t$end\t$ref\t$alt\t$gene\t$within\t$donGain\t$accGain\t$donLoss\t$accLoss")
#        fi
#    done < output/"$fileID"_out.txt
#    echo -e "$prevLine" >> output/"$fileID"_out.txt_merged
#end comment to slip duplicate operations
done
