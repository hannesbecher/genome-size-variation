# requires kmc in path

# -ci1 is important because this is an assembly. Most k-mer would be discarded otherwise.
kmc -k21 -fm -t10 -ci1 Option_1_A1.fasta plastidDB21 .
kmc -k27 -fm -t10 -ci1 Option_1_A1.fasta plastidDB27 .
