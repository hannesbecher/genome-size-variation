
# script to generate k-mer DBs and "full-range" spectra
#script requires KMC3 and kmc_tools installed

WORK_DIR=$(pwd)

# replace with your own sample names
samples=(E030 E031 E032 E040 E065 E068 E073)

for sample in "${samples[@]}"
do
	echo "Currently creating histograms for $sample"

	# path to where the trimmed read files are, ADJUST!
	ls /disk2/hbecher_tmp/trimmedFastp/$sample*gz > file_names

	kmc -t20 -k21 -cs150000000 @file_names $sample"full" . #generate kmer dataset, high multiplicity, requires file_names of input genomes
	rm file_names
	kmc_tools transform $sample"full" histogram $sample"full.hist" -cx150000000 #generate histogram table from kmer dataset
	awk '{ if( $2 != 0 ){ print $0 } }' $sample"full.hist" > $sample"full.hist.no0" #remove lines with zeros from histogram
	rm $sample"full.hist"

	# path to plastid database, ADJUST!
	kmc_tools simple $sample"full" /disk2/hbecher_tmp/GSproject/plastid_genome/plastidDB21 kmers_subtract $sample"noPl" -cs150000000 #generate new kmer dataset with kmers corresponding to plastid genome removed

	rm $sample"full.kmc*"
	kmc_tools transform $sample"noPl" histogram $sample"noPl.hist" -cx150000000
	awk '{if ($2 != 0){print $0}}' $sample"noPl.hist" > $sample"noPl.hist.no0"
	rm $sample"noPl.hist"

	# path to mito database, ADJUST!
	kmc_tools simple $sample"noPl" /disk2/hbecher_tmp/GSproject/mito_genome/mitoDB21 kmers_subtract $sample"noPlMt" -cs150000000
	rm $sample"noPl.kmc*"2
	kmc_tools transform $sample"noPlMt" histogram $sample"noPlMt.hist" -cx150000000 #generate new histogram with kmers corresponding to plastid genome removed
	#rm $sample"noPlMt.kmc*"
	awk '{ if( $2 != 0 ){ print $0 } }' $sample"noPlMt.hist" > $sample"noPlMt.hist.no0" #generate new histogram with no plastid kmers and no zeros
	rm $sample"noPlMt.hist"
	echo "Histograms for $sample completed"

done
