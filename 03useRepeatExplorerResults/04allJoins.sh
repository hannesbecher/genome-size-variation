# carry out dumping and joining for all diploids with reference sample E030

samples=("E031" "E032" "E040" "E065" "E068" "E073")
ndeps=(42.4 35 67.4 28.5 25.5 20.8) # monoploid depths from Tetmer (after fastp trimming)
d30=54 # monoploid depth of refeecne sample E030

# loop over sample indices
for i in {0..5}; do
    echo ${samples[i]}
    echo "Dumping..."
    # loop over super clusters
    for j in {01..50}; do
        kmc_tools transform ${samples[i]}/intersect${samples[i]}RE$j dump -s ${samples[i]}/${samples[i]}RE$j.dump

        echo "Joining..."
        join -a 1 -a 2 -e 0 -o 0,1.2,2.2 E030/E030RE$j.dump ${samples[i]}/${samples[i]}RE$j.dump >  ${samples[i]}/${samples[i]}RE$j.join && rm ${samples[i]}/${samples[i]}RE$j.dump

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
        ' ${samples[i]}/${samples[i]}RE$j.join | sort | uniq -c > ${samples[i]}/${samples[i]}RE$j.binCounts && rm  ${samples[i]}/${samples[i]}RE$j.join

        echo "Removing leading whitespace..."
        sed "s/^[ \t]*//" -i ${samples[i]}/${samples[i]}RE$j.binCounts

        echo "Done."
    done
done

echo "ALL DONE."

