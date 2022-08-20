#!/usr/bin/env bash

# Version 2.8 for use with large genome resequencing projects, specifically the eastern oyster

export LC_ALL=en_US.UTF-8
export SHELL=bash

#mtDNA Chromosome
mtDNA="NC_007175.2"

# This script uses the individual raw.bcf files for filtering instead of a single VCF file.  This enables parallelization and much faster processing of larger files

# Usage is `PVCF_filter PREFIX PROCESSORS`

if [[ -z "$2" ]]; then
echo "Usage is bash PVCF_filter PREFIX PROCESSORS"
exit 1
fi

if [[ -z "$3" ]]; then
PE="no"
else
PE=$3
fi

NumProc=$2
PREFIX=$1
SProc=$(($NumProc / 5))

NumInd=$(ls raw.*.vcf.gz 2> /dev/null | wc -l)
NumInd=$(($NumInd - 0))





if [ "$NumInd" -gt 0 ];then
	GZIP="TRUE"
fi

if [ "$GZIP" == "TRUE" ]; then
	ls raw.*.vcf.gz | sed 's/raw.//g' | sed 's/.vcf.gz//g' > list
	cat list | parallel --no-notice -j $NumProc "bcftools view -i 'F_MISSING<0.75 && QUAL > 29' raw.{}.vcf.gz -O v | bcftools +setGT -- -t q -n . -i 'FORMAT/DP<5' 2> $PREFIX.filter.errors | bcftools view -i 'F_MISSING<0.25 && MAF > 0.0001' -O b -o $PREFIX.TRSdp5g75.{}.recode.bcf" 2>> $PREFIX.filter.errors
	FIRST=$(head -1 list)
	bcftools view raw.$FIRST.vcf.gz | head -1000 | grep "##contig" | sed 's/##contig=<ID=//g' | sed 's/,length=/	1	/g' | sed 's/>//g' > $PREFIX.genome.file
else
	NumInd=$(ls raw.*.bcf 2> /dev/null | wc -l)
	NumInd=$(($NumInd - 0))
	ls raw.*.bcf | sed 's/raw.//g' | sed 's/.bcf//g' > list
	cat list | parallel --no-notice -j $NumProc "bcftools view -i 'F_MISSING<0.75 && QUAL > 29' raw.{}.bcf  | bcftools +setGT -- -t q -n . -i 'FORMAT/DP<5' 2> $PREFIX.filter.errors | bcftools view -i 'F_MISSING<0.25 && MAF > 0.0001' -O b -o $PREFIX.TRSdp5g75.{}.recode.bcf 2>> $PREFIX.filter.errors"
	FIRST=$(head -1 list)
	bcftools view raw.$FIRST.bcf | head -1000 | grep "##contig" | sed 's/##contig=<ID=//g' | sed 's/,length=/	1	/g' | sed 's/>//g' > $PREFIX.genome.file

fi

NumK=$(($NumInd/1000))

if [ ! -d "raw" ]; then
	mkdir raw
fi

mv raw.*.vcf.gz raw.*.bcf ./raw 2> /dev/null

empty_test(){
	PRE=$2
	POST=$3
	SNPs=$(bcftools query -f '%CHROM  %POS\n' $PRE.$1.$POST | wc -l)
        if [ "$SNPs" -lt 1 ];then
        	rm $PRE.$1.$POST
	else
		echo $1
        fi
	}
	
export -f empty_test

cat list | parallel --env empty_test --no-notice -j $NumProc -k "empty_test {} $PREFIX.TRSdp5g75 recode.bcf" > list.no.empty

mv list.no.empty list

cat list | parallel --no-notice -j $NumProc "bcftools index $PREFIX.TRSdp5g75.{}.recode.bcf 2>> $PREFIX.filter.errors"

#if [ ! -d "TRSdp10" ]; then
#	mkdir TRSdp10
#fi

#mv $PREFIX.TRSdp5.*.recode.vcf.gz TRSdp10

#mtDNA processing
FIRST=$(head -1 list)

grep $mtDNA $PREFIX.genome.file > $PREFIX.mtDNA.file
grep -v $mtDNA $PREFIX.genome.file > $PREFIX.nDNA.file


cat list | parallel --no-notice -j $NumProc "bcftools view -R $PREFIX.nDNA.file $PREFIX.TRSdp5g75.{}.recode.bcf -O b -o $PREFIX.TRSdp5g75.nDNA.{}.bcf"
cat list | parallel --no-notice -j $NumProc "bcftools view -R $PREFIX.mtDNA.file $PREFIX.TRSdp5g75.{}.recode.bcf -O v -o $PREFIX.TRSdp5g75.mtDNA.{}.vcf"

cat list | parallel --env empty_test --no-notice -j $NumProc -k "empty_test {} $PREFIX.TRSdp5g75.mtDNA vcf" > list.mtdna
	
cat list.mtdna | parallel --no-notice -j $NumProc "vcffilter -s -f 'AB < 0.001' $PREFIX.TRSdp5g75.mtDNA.{}.vcf | vcffilter -s -f 'QUAL / DP > 0.25' > $PREFIX.TRSdp5g75.mtDNA.F.{}.vcf "
cat list.mtdna | parallel --no-notice -j $NumProc "vcfallelicprimitives -k -g $PREFIX.TRSdp5g75.mtDNA.F.{}.vcf | sed 's:\.|\.:\.\/\.:g' | bcftools view --threads 2 -V indels,other -O v -o $PREFIX.SNP.TRSdp5g75.mtDNA.{}.vcf 2>> $PREFIX.filter.errors"
																																
vcfcombine $PREFIX.SNP.TRSdp5g75.mtDNA.*.vcf | vcfstreamsort -a > $PREFIX.SNP.TRSdp5g75.mtDNA.vcf

cat list.mtdna | parallel --no-notice -j $NumProc "bgzip $PREFIX.SNP.TRSdp5g75.mtDNA.{}.vcf"
cat list.mtdna | parallel --no-notice -j $NumProc "bgzip $PREFIX.TRSdp5g75.mtDNA.F.{}.vcf"
cat list.mtdna | parallel --no-notice -j $NumProc "bgzip $PREFIX.TRSdp5g75.mtDNA.{}.vcf"


if [ ! -d "mtDNA.vcf" ]; then
	mkdir mtDNA.vcf
fi

if [ ! -d "filtered" ]; then
	mkdir filtered
fi

mv $PREFIX.SNP.TRSdp5g75.mtDNA.*.vcf.gz $PREFIX.TRSdp5g75.mtDNA.*.vcf.gz ./mtDNA.vcf

mv $PREFIX.SNP.TRSdp5g75.mtDNA.vcf ./filtered


#if [ ! -d "TRSdp5g75" ]; then
#        mkdir TRSdp5g75
#fi

#mv $PREFIX.TRSdp5g75.*.recode.vcf.gz $PREFIX.TRSdp5g75.*.recode.vcf.gz.tbi TRSdp5g75
 
rm $PREFIX.TRSdp5g75.*.recode.bcf*
  
#nDNA processing
echo "This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,"
echo "quality versus depth, allelic balance at heterzygous individuals, and paired read representation."
echo -e "The script assumes that loci and individuals with low call rates (or depth) have already been removed. \n"
echo -e "Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters \n"

PREFIX="$PREFIX.TRSdp5g75.nDNA"

#Creates a file with the original site depth and qual for each locus

AWK1='!/#/ {print $1 "\t" $2 "\t" $6}'
AWK2='!/#/ {print $1}'
AWK3='!/NP/ && !/#/'
AWK4='!/#/ {print $1 "\t" $2}'

cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.bcf | cut -f8  | grep -P -oe 'DP=[0-9]*' | sed -s 's/DP=//g' > TEMP.{}.DEPTH"

cat TEMP.*.DEPTH > $PREFIX.DEPTH
rm TEMP*.DEPTH

cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.bcf | mawk '$AWK1' " > TEMP.{}.loci.qual
cat TEMP.*.loci.qual > $PREFIX.loci.qual
rm TEMP.*.loci.qual

#Filters out sites with that on average have heterzygotes with less than a 0.05 allele balance between reads from each allele and Quality / Depth < 0.1
cat list | parallel --no-notice -j $NumProc "tabix $PREFIX.{}.bcf"
cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.bcf | vcffilter -s -f 'AB > 0.1 & AB < 0.9 | AB < 0.01 | AB > 0.99'  | vcffilter -s -f 'AF < 0.01 | QUAL / DP > 0.1' | bcftools view -O b -o $PREFIX.{}.temp.bcf "

FILTERED=$( cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
OLD=$( cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )

NUMFIL=$(($OLD - $FILTERED))
echo -e "Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth\n" $NUMFIL "of" $OLD "\n"
echo -e "Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth\n" $NUMFIL "of" $OLD "\n" > $PREFIX.filterstats

if [ "$PE" != "yes" ]; then
	FILTERED2=$(cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | mawk '$AWK2' | wc -l " | mawk '{sum = sum + $1} END {print sum}' )
	NUMFIL2=$(($FILTERED - $FILTERED2))
	echo -e "Number of additional sites filtered based on properly paired status\n" $NUMFIL2 "of" $FILTERED "\n"
	echo -e "Number of additional sites filtered based on properly paired status\n" $NUMFIL2 "of" $FILTERED "\n" >> $PREFIX.filterstats
else	
	FILTERED2=$(cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A | mawk '$AWK3' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
	NUMFIL2=$(($FILTERED - $FILTERED2))
	echo -e "Number of additional sites filtered based on properly paired status\n" $NUMFIL2 "of" $FILTERED "\n"
	echo -e "Number of additional sites filtered based on properly paired status\n" $NUMFIL2 "of" $FILTERED "\n" >> $PREFIX.filterstats

fi

wait
wait
wait

#Calculates the average depth and three times the square root of mean depth
DEPTH=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $PREFIX.DEPTH)
MD=$(python -c "import math; print( round(math.sqrt($DEPTH)))")
MD=$(python -c "import math; print( round( $MD * 3))")
DEPTH=$(python -c "import math; print( round($DEPTH + $MD))")

#Filters loci above the mean depth + 1 standard deviation that have quality scores that are less than 2*DEPTH
paste $PREFIX.loci.qual $PREFIX.DEPTH | mawk -v x=$DEPTH '$4 > x'| mawk -v x=$DEPTH '$3 < 2 * x' > $PREFIX.lowQDloc
cut -f1,2 $PREFIX.lowQDloc > $PREFIX.lowQDloci

SITES=$(cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcftools --vcf - --exclude-positions $PREFIX.lowQDloci 2>&1 | grep Sites"| cut -f4 -d " " | mawk '{sum = sum + $1} END {print sum}' )
LQDL=$(( $FILTERED - $SITES ))

echo -e "Number of sites filtered based on high depth and lower than 2*DEPTH quality score\n" $LQDL "of" $FILTERED "\n"
echo -e "Number of sites filtered based on high depth and lower than 2*DEPTH quality score\n" $LQDL "of" $FILTERED "\n" >> $PREFIX.filterstats

#Recalculates site depth for sites that have not been previously filtered
NPFILTER='FILTER!="NP"'
DQUERY='%INFO/DP\n'
if [ "$PE" != "yes" ]; then

	cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | bcftools view -f .,PASS -T ^$PREFIX.lowQDloci | bcftools query -f '$DQUERY' > $PREFIX.{}.ldepth 2>> $PREFIX.filter.errors"
else
	cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A  | bcftools view -f .,PASS -T ^$PREFIX.lowQDloci | bcftools query -f '$DQUERY' > $PREFIX.{}.ldepth 2>> $PREFIX.filter.errors"
fi

cat list | parallel --no-notice -k -j $NumProc "cut -f1 $PREFIX.{}.ldepth" | mawk '!/SUM_DEPTH/' > TEMP.{}.site.depth
cat TEMP.*.site.depth > $PREFIX.site.depth
rm TEMP.*.site.depth

cat list | parallel --no-notice -k -j $NumProc gzip $PREFIX.{}.ldepth

DP=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $PREFIX.site.depth)
SD=$(mawk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $PREFIX.site.depth)

#Calculates actual number of individuals in VCF file
#This is important because loci will now be filtered by mean depth calculated with individuals present in VCF
FFILE=$(head -1 list)
IND=$(bcftools view $PREFIX.$FFILE.temp.bcf |head -1000 | mawk '/#/' | tail -1 | wc -w)
IND=$(($IND - 9))

mawk '!/D/' $PREFIX.site.depth | mawk -v x=$IND '{print $1/x}' > meandepthpersite

#Calculates a mean depth cutoff to use for filtering
DP=$(perl -e "print ($DP+ 1.645*$SD) / $IND")
PP=$(mawk '!/SUM/' $PREFIX.site.depth | sort -rn | perl -e '$d=.05;@l=<>;print $l[int($d*$#l)]' )
PP=$(perl -e "print int($PP / $IND)")
GP=$(perl -e "print int($PP * 1.25)")
export GP

gnuplot << \EOF >> $PREFIX.filterstats
set terminal dumb size 120, 30
set autoscale
high=system("echo $GP")
set xrange [10:high]
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
#set yr [0:100000]
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics floor(high/20)
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
high=system("echo $GP")
set xrange [10:high]
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
#set yr [0:100000]
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics floor(high/20)
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


if [[ -z "$4" ]]; then
echo "The 95% cutoff would be" $PP
echo "Would you like to use a different maximum mean depth cutoff than "$PP", yes or no"

#read NEWCUTOFF
NEWCUTOFF=no
else
NEWCUTOFF=$4
fi

if [ "$NEWCUTOFF" != "yes" ]; then
echo -e "Maximum mean depth cutoff is" $PP
echo -e "Maximum mean depth cutoff is" $PP >> $PREFIX.filterstats

else
	if [[ -z "$5" ]]; then
		echo "Please enter new cutoff"
		read PP
	else
		PP=$5
	fi
echo -e "Maximum mean depth cutoff is" $PP >> $PREFIX.filterstats
fi

#Combines all filters to create filtered VCF files
if [ "$PE" != "yes" ]; then

	cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | vcftools --vcf - --remove-filtered NP --max-meanDP $PP --recode --exclude-positions $PREFIX.lowQDloci --recode-INFO-all --stdout 2>> /dev/null | vcfstreamsort -a | bcftools view -O b -o $PREFIX.{}.FIL.bcf 2>> $PREFIX.filter.errors"
else
	cat list | parallel --no-notice -k -j $NumProc "bcftools view $PREFIX.{}.temp.bcf | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A  | vcftools --vcf - --remove-filtered NP --max-meanDP $PP --recode --exclude-positions $PREFIX.lowQDloci --recode-INFO-all --stdout 2>> /dev/null | vcfstreamsort -a | bcftools view -O b -o $PREFIX.{}.FIL.bcf 2>> $PREFIX.filter.errors"

fi

rm $PREFIX.*.temp.bcf

FILTERED3=$(cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.FIL.bcf -O v | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' ) && OLD2=$(cat $PREFIX.site.depth | wc -l)

NUMFIL3=$(($OLD2 - $FILTERED3))

echo -e "Number of sites filtered based on maximum mean depth\n" $NUMFIL3 "of" $OLD2 "\n"
echo -e "Number of sites filtered based on maximum mean depth\n" $NUMFIL3 "of" $OLD2 "\n" >> $PREFIX.filterstats

NUMFIL4=$(($OLD - $FILTERED3))

echo -e "Total number of sites filtered\n" $NUMFIL4 "of" $OLD "\n"
echo -e "Total number of sites filtered\n" $NUMFIL4 "of" $OLD "\n" >> $PREFIX.filterstats

echo -e "Remaining sites\n" $FILTERED3 "\n"
echo -e "Remaining sites\n" $FILTERED3 "\n" >> $PREFIX.filterstats


if [ ! -d "nDNA" ]; then
	mkdir nDNA
	mkdir metrics
	mkdir SNPs
	mkdir nDNA.INDels
#	mkdir nDNA.FIL
fi

mv $PREFIX.[0-9]*[!FIL].bcf $PREFIX.[0-9]*.bcf.csi ./nDNA


cat list | parallel --env empty_test --no-notice -j $NumProc -k "empty_test {} $PREFIX FIL.bcf" > list.no.empty

mv list.no.empty list

echo -e "Variants will now be composed into SNPs and INDels\n"
echo -e "Variants will now be composed into SNPs and INDels\n" >> $PREFIX.filterstats

#cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.FIL.bcf | vcfallelicprimitives -k -g | sed 's:\.|\.:\.\/\.:g' | bgzip -c > $PREFIX.{}.prim.vcf.gz"

cat list | parallel --no-notice -j $NumProc "bcftools view $PREFIX.{}.FIL.bcf | vcfallelicprimitives -k -g | sed 's:\.|\.:\.\/\.:g' | bcftools +fill-tags -O b -o $PREFIX.{}.prim.bcf"

#mv $PREFIX.*.FIL.bcf ./nDNA.FIL/ 2>> $PREFIX.filter.errors

rm $PREFIX.*.FIL.bcf

cat list | parallel --no-notice -j $NumProc "bcftools view --threads 2 -v indels,other $PREFIX.{}.prim.bcf -O z -o INDELS.$PREFIX.{}.bcf" 2>> $PREFIX.filter.errors 
cat list | parallel --no-notice -j $NumProc "bcftools view --threads 2 -V indels,other $PREFIX.{}.prim.bcf -O z -o SNP.$PREFIX.{}.FIL.bcf" 2>> $PREFIX.filter.errors 

rm $PREFIX.*.prim.bcf

FILTERED4=$(cat list | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
FILTERED5=$(cat list | parallel --no-notice -j $NumProc "bcftools view INDELS.$PREFIX.{}.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )

echo -e "Number of SNPs before filtering for call rate and minor allele frequency\n" $FILTERED4 "\n"
echo -e "Number of SNPs before filtering for call rate and minor allele frequency\n" $FILTERED4 "\n" >> $PREFIX.filterstats

echo -e "Number of INDels\n" $FILTERED5 "\n"
echo -e "Number of INDels\n" $FILTERED5 "\n" >> $PREFIX.filterstats

mv INDELS.$PREFIX.*.bcf ./nDNA.INDels 2>> $PREFIX.filter.errors

# Filter by genotype and minor allele frequency

empty_loci(){
	PRE=$2
	SNPs=$(cat $PRE.$1 | wc -l)
        if [ "$SNPs" -lt 1 ];then
        	rm $PRE.$1
	else
		echo $1
        fi
	}
export -f empty_loci

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags |bcftools view -i 'F_MISSING < 0.1 & MAF > 0.01 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g9.maf01.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED6=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g9.maf01.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 90% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED6 "\n"
echo -e "Number of SNPs with 90% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED6 "\n" >> $PREFIX.filterstats


cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -i 'F_MISSING < 0.1 & MAF > 0.05 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g9.maf05.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED7=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g9.maf05.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 90% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED7 "\n"
echo -e "Number of SNPs with 90% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED7 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -i 'F_MISSING < 0.05 & MAF > 0.01 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g95.maf01.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED8=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g95.maf01.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 95% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED8 "\n"
echo -e "Number of SNPs with 95% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED8 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -i 'F_MISSING < 0.05 & MAF > 0.05 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g95.maf05.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED9=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g95.maf05.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 95% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED9 "\n"
echo -e "Number of SNPs with 95% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED9 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -i 'F_MISSING = 0 & MAF > 0.01 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g1.maf01.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED10=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g1.maf01.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 100% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED10 "\n"
echo -e "Number of SNPs with 100% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED10 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -i 'F_MISSING = 0 & MAF > 0.05 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g1.maf05.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED11=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g1.maf05.{}.FIL.bcf | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of SNPs with 100% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED11 "\n"
echo -e "Number of SNPs with 100% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED11 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -m 2 -M 2 -i 'F_MISSING = 0 & MAF > 0.01 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g1.maf01.max2alleles.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED13=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g1.maf01.max2alleles.{}.FIL.bcf 2>/dev/null | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of biallelic SNPs with 100% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED13 "\n"
echo -e "Number of biallelic SNPs with 100% call rate, MAF of 0.01, and excess heterozygous loci removed\n" $FILTERED13 "\n" >> $PREFIX.filterstats

cat list | parallel --no-notice -j $NumProc "bcftools view -S filter.set SNP.$PREFIX.{}.FIL.bcf 2>/dev/null | bcftools +fill-tags  |bcftools view -m 2 -M 2 -i 'F_MISSING = 0 & MAF > 0.05 & ExcHet > 1.0e-8' -O v 2>/dev/null | mawk '$AWK4' > loci.{}" 2>> $PREFIX.filter.errors
cat list | parallel --env empty_loci --no-notice -j $NumProc -k "empty_loci {} loci " > list.loci
cat list.loci | parallel --no-notice -j $NumProc "bcftools view -T loci.{} SNP.$PREFIX.{}.FIL.bcf -O b -o SNP.$PREFIX.g1.maf05.max2alleles.{}.FIL.bcf" 2>> $PREFIX.filter.errors
FILTERED12=$(cat list.loci | parallel --no-notice -j $NumProc "bcftools view SNP.$PREFIX.g1.maf05.max2alleles.{}.FIL.bcf 2>/dev/null | mawk '$AWK2' | wc -l" | mawk '{sum = sum + $1} END {print sum}' )
cat list.loci | parallel --no-notice -j $NumProc rm loci.{}
echo -e "Number of biallelic SNPs with 100% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED12 "\n"
echo -e "Number of biallelic SNPs with 100% call rate, MAF of 0.05, and excess heterozygous loci removed\n" $FILTERED12 "\n" >> $PREFIX.filterstats

# Collate
seq 0 $NumK | parallel --no-notice -j $NumProc "ls SNP.$PREFIX.0{}*.FIL.bcf > bcf.{}.list"
seq 0 $NumK | parallel --no-notice -j $NumProc "bcftools concat -n -f bcf.{}.list -O b -o SNP.$PREFIX.Collated.{}.FIL.bcf 2>> $PREFIX.filter.errors"

ls SNP.$PREFIX.Collated.*.FIL.bcf > collated.bcf.list
bcftools concat -n -f collated.bcf.list -O b 2>> $PREFIX.filter.errors | bcftools view -O z --threads $NumProc -o SNP.$PREFIX.FIL.vcf.gz

rm SNP.$PREFIX.Collated.*.FIL.bcf

mv SNP.$PREFIX.[0-9]*.FIL.bcf SNP.$PREFIX.FIL.vcf.gz ./SNPs

for i in {g9,g95,g1};
do
	for j in {maf01,maf05};
	do
	seq 0 $NumK | parallel --no-notice -j $NumProc "ls SNP.$PREFIX.$i.$j.0{}*.FIL.bcf > bcf.{}.list"
	seq 0 $NumK | parallel --no-notice -j $SProc "bcftools concat -n -f bcf.{}.list -O b -o SNP.$PREFIX.$i.$j.Collated.{}.FIL.bcf 2>> $PREFIX.filter.errors"
	rm SNP.$PREFIX.$i.$j.[0-9]*.FIL.bcf
	ls SNP.$PREFIX.$i.$j.Collated.*.FIL.bcf > collated.bcf.list
	bcftools concat -n -f collated.bcf.list -O b 2>> $PREFIX.filter.errors | bcftools view -O z --threads $NumProc -o SNP.$PREFIX.$i.$j.FIL.vcf.gz
	rm SNP.$PREFIX.$i.$j.Collated.*.FIL.bcf
	done
done
	
seq 0 $NumK | parallel --no-notice -j $NumProc "ls SNP.$PREFIX.g1.maf05.max2alleles.0{}*.FIL.bcf > bcf.{}.list"
seq 0 $NumK | parallel --no-notice -j $SProc "bcftools concat -n -f bcf.{}.list -O b -o SNP.$PREFIX.Collated.g1.maf05.max2alleles.{}.FIL.bcf 2>> $PREFIX.filter.errors"

rm SNP.$PREFIX.g1.maf05.max2alleles.[0-9]*.FIL.bcf
ls SNP.$PREFIX.Collated.g1.maf05.max2alleles.*.FIL.bcf > collated.bcf.list
bcftools concat -n -f collated.bcf.list -O b 2>> $PREFIX.filter.errors | bcftools view -O z --threads $NumProc -o SNP.$PREFIX.g1.maf05.max2alleles.FIL.vcf.gz

rm bcf.*.list
rm collated.bcf.list
rm SNP.$PREFIX.Collated.g1.maf05.max2alleles.*.FIL.bcf

seq 0 $NumK | parallel --no-notice -j $NumProc "ls SNP.$PREFIX.g1.maf01.max2alleles.0{}*.FIL.bcf > bcf.{}.list"
seq 0 $NumK | parallel --no-notice -j $SProc "bcftools concat -n -f bcf.{}.list -O b -o SNP.$PREFIX.Collated.g1.maf01.max2alleles.{}.FIL.bcf 2>> $PREFIX.filter.errors"

rm SNP.$PREFIX.g1.maf01.max2alleles.[0-9]*.FIL.bcf
ls SNP.$PREFIX.Collated.g1.maf01.max2alleles.*.FIL.bcf > collated.bcf.list
bcftools concat -n -f collated.bcf.list -O b 2>> $PREFIX.filter.errors | bcftools view -O z --threads $NumProc -o SNP.$PREFIX.g1.maf01.max2alleles.FIL.vcf.gz
	
rm bcf.*.list
rm collated.bcf.list
rm SNP.$PREFIX.Collated.g1.maf01.max2alleles.*.FIL.bcf

mv SNP.$PREFIX.*.*.FIL.vcf.gz ./filtered

ls $PREFIX.DEPTH $PREFIX.lo* $PREFIX.site.depth meandepthpersite | parallel --no-notice -j $NumProc gzip {}
mv $PREFIX.DEPTH.gz $PREFIX.lo* $PREFIX.site.depth.gz meandepthpersite.gz  ./metrics
cat $PREFIX.*.ldepth.gz > $PREFIX.ldepth.gz 
mv $PREFIX.ldepth.gz ./metrics
rm $PREFIX.*.ldepth.gz

mv $PREFIX.filterstats ./filtered

echo -e "Filter stats and filtered VCF files stored in $PREFIX.filterstats\n"
echo -e "Both are stored in the filtered directory\n"

