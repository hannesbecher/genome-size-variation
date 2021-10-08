# Script to bin multiplicities in joined k-mer dumps
# Assumes three columns, applies the function normBin to columns 2 and 3.
# USAGE bash scale_dum.sh <joined dump> <noramalisation factor 1> <normalisation factor 2>

# The normalisation factors for spectra are read from the command line:
fact1=$2
fact2=$3
# NB $1, $2, and $3 have a different meaning inside the awk string below.

# Location of the dump file
dfile=$1
#echo $dfile

# A binning function
awk -v f1=$fact1 -v f2=$fact2 '
function normBin(num, factor){
    if(num/factor < 0.5)
    return 0 
    return int(log(num/factor/0.5)/log(1.1))
}
{
    print normBin($2, f1) " " normBin($3, f2)
}
' $dfile
