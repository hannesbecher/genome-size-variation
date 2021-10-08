# Generate bin counts from join
echo "Transforming dump..."
bash scale_dump.sh <(zcat testJoin.gz) 15 13 | sort | uniq -c > transformedJoin

# Remove leading whitespace
echo "Removing leading white space..."
sed "s/^[ \t]*//" -i transformedJoin

# Plot histograms
echo "Plotting histograms..."
Rscript hist2D.R transformedJoin

echo "Done."
