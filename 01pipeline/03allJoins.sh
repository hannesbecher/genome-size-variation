# carry out dumping and joining for all diploids with reference sample E030

samples=("E031" "E032" "E040" "E065" "E068" "E073")
ndeps=(42.4 35 67.4 28.5 25.5 20.8) # depths after fastp trimming!
d30=54 # post trimming

#kmc_tools transform E030noPlMt dump -s E030noPlMt.dump

for i in {0..5}; do
    echo ${samples[i]}
    echo "Dumping..."
    kmc_tools transform ${samples[i]}noPlMt dump -s ${samples[i]}noPlMt.dump

    echo "Joining..."
    join -a 1 -a 2 -e 0 -o 0,1.2,2.2 E030noPlMt.dump ${samples[i]}noPlMt.dump > E030${samples[i]}noPlMt.join && rm ${samples[i]}noPlMt.dump

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
    ' E030${samples[i]}noPlMt.join | sort | uniq -c > E030${samples[i]}noPlMt.binCounts.trim && rm E030${samples[i]}noPlMt.join

    echo "Removing leading whitespace..."
    sed "s/^[ \t]*//" -i E030${samples[i]}noPlMt.binCounts.trim

    echo "Done."
done

echo "ALL DONE."

