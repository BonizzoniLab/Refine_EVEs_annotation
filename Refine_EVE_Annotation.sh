#!/bin/bash -x

#----------------------------------------------------------------------------------------------------------------------------
# Input settings
#----------------------------------------------------------------------------------------------------------------------------
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -pipeline_directory)
    PIPELINE_DIR="$2"
    shift # past argument
    ;;
    -reverse_blast_tool)
    TOOL="$2"
    shift # past argument
    ;;
	-file_blastx)
	BLASTX="$2"
	shift # past argument
	;;
	-file_bed_tophit)
	TOPHIT="$2"
	shift # past argument
	;;
    -VHC_directory)
	VHCDIR="$2"
	shift # past argument
	;;
    -output_directory)
	OUTDIR="$2"
	shift # past argument
	;;
    -taxonkit_exe)
	TAXONKIT_EXE="$2"
	shift # past argument
	;;
	--default)
	DEFAULT=YES
	;;
	*)
	# unknown option
	;;			
esac
shift # past argument or value
done

START=$(date +%s);

#----------------------------------------------------------------------------------------------------------------------------
# Control settings
#----------------------------------------------------------------------------------------------------------------------------
if [ -z "$TOOL" ];then
  echo "Please select the tool used for the reverse blastx: 'blastx' or 'diamond'. Default: 'blastx'."
  TOOL="blastx"
fi

#----------------------------------------------------------------------------------------------------------------------------
# Creation of the directories
#----------------------------------------------------------------------------------------------------------------------------

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
    mkdir $OUTDIR/Output
	mkdir $OUTDIR/VHC
else
    if [ ! -d $OUTDIR/Output ]; then
        mkdir $OUTDIR/Output
    fi
    if [ ! -d $OUTDIR/VHC ]; then
        mkdir $OUTDIR/VHC
    fi
fi
printf "\n"
echo -e $OUTDIR 'Created'

#----------------------------------------------------------------------------------------------------------------------------
# Extraction of taxid from the blastx or diamond file to classified them as viral or non-viral using VHC
#----------------------------------------------------------------------------------------------------------------------------

if [ "$TOOL" = "diamond" ]; then
	echo -e 'The tool used for the reverse blastx is diamond'
    python $PIPELINE_DIR/SelectUniqueTaxid.py \
    -i_Blastx $BLASTX \
    -position 11 \
    > $VHCDIR/TaxonID.tsv
	
else
	echo -e 'The tool used for the reverse blastx is blastx'
    python $PIPELINE_DIR/SelectUniqueTaxid.py \
    -i_Blastx $BLASTX \
    -position 12 \
    > $VHCDIR/TaxonID.tsv

fi

#----------------------------------------------------------------------------------------------------------------------------
# Classification of taxid as viral or non-viral using VHC
#----------------------------------------------------------------------------------------------------------------------------

cd $VHCDIR

python3 VHost_Classifier.py TaxonID.tsv virushostdb.tsv $OUTDIR/VHC

cd $OUTDIR

#----------------------------------------------------------------------------------------------------------------------------
# Assign to each tophit the best viral and non viral result from blastx
#----------------------------------------------------------------------------------------------------------------------------

if [ "$TOOL" = "diamond" ]; then

    python $PIPELINE_DIR/CheckVHclassifier_diamond.py \
    -i_VHname $OUTDIR/VHC/Virus.csv \
    -i_Blastx $BLASTX \
    -i_TopHits $TOPHIT \
    -output_path $OUTDIR/Output
	
else
    python $PIPELINE_DIR/CheckVHclassifier.py \
    -i_VHname $OUTDIR/VHC/Virus.csv \
    -i_Blastx $BLASTX \
    -i_TopHits $TOPHIT \
    -output_path $OUTDIR/Output

fi

#----------------------------------------------------------------------------------------------------------------------------
# Extraction of the taxid from the best results
#----------------------------------------------------------------------------------------------------------------------------

python $PIPELINE_DIR/SelectUniqueTaxid.py \
-i_Blastx $OUTDIR/Output/CompleteTable.txt \
-position 19 \
> $OUTDIR/Output/Viral_taxid.txt

#----------------------------------------------------------------------------------------------------------------------------
# Extraction of the viral family and the viral order of the selected taxid
#----------------------------------------------------------------------------------------------------------------------------

$TAXONKIT_EXE lineage -t $OUTDIR/Output/Viral_taxid.txt > $OUTDIR/Output/Taxonkit_lineage_Viral_taxid.txt

$TAXONKIT_EXE reformat -t -f "{o};{f}" $OUTDIR/Output/Taxonkit_lineage_Viral_taxid.txt > $OUTDIR/Output/Taxonkit_reformat_Order-Family_Viral_taxid.txt

#----------------------------------------------------------------------------------------------------------------------------
# Generate the output file
#----------------------------------------------------------------------------------------------------------------------------

python $PIPELINE_DIR/AssignOrderFamily_viralHits.py \
-CompleteTable $OUTDIR/Output/CompleteTable.txt \
-Taxonkit_out $OUTDIR/Output/Taxonkit_reformat_Order-Family_Viral_taxid.txt


printf "\n"     		
echo -e 'Time required'
END=$(date +%s);
echo $((END-START)) | awk '{print int($1/60)":"int($1%60)" min:sec"}'
echo -e 'End'