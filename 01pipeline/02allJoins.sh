# carry out dumping and joining for all diploids with reference sample E030
ref="E030" # reference individual
samples=("E031" "E032" "E040" "E065" "E068" "E073") # other individuals to analyse
ndeps=(42.4 35 67.4 28.5 25.5 20.8) # per-individual depths of the heterozygous peak (after fastp trimming)
d30=54 # reference depths of the heterozygous peak

kmc_tools transform $ref'noPlMt' dump -s $ref'noPlMt.dump'

for i in {0..5}; do
    echo ${samples[i]}
    echo "Dumping..."
    kmc_tools transform ${samples[i]}noPlMt dump -s ${samples[i]}noPlMt.dump

    echo "Joining..."
    join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $ref'noPlMt.dump' ${samples[i]}noPlMt.dump > $ref${samples[i]}noPlMt.join && rm ${samples[i]}noPlMt.dump

    echo "Scaling and binning..."
    awk -v f1=$d30 -v f2=${ndeps[i]} '
    function normBin(num, factor){
    if(num/factor < 0.5)
    return 0
    return int(log(num/factor/0.5)/log(1.1)) + 1
    }
    {
    print normBin($2, f1) " " normBin($3, f2)
    }
    ' $ref${samples[i]}noPlMt.join | sort | uniq -c > $ref${samples[i]}noPlMt.binCounts.trim && rm $ref${samples[i]}noPlMt.join

    echo "Removing leading whitespace..."
    sed "s/^[ \t]*//" -i $ref${samples[i]}noPlMt.binCounts.trim

    echo "Done."
done
rm $ref'noPlMt.dump'
echo "ALL DONE."

