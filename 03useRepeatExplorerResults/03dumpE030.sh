# make sure kmc is in path
# This is to dump the k-mer or the reference sample to which each sample is going to be compared.
# The other sample k-mers are dumped only temporarily with the next script.
for i in {01..50}; do kmc_tools transform E030/intersectE030RE$i dump -s E030/E030RE$i.dump; done
