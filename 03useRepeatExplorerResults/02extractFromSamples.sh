# make sure kmc_tools is in path

pref=""

# loop over sample names (replace with your sample names)
for s in E030 E031 E032 E040 E065 E068 E073; do
  echo $s
  mkdir $s
  # loop over suple cluster numbers
  for i in {01..50}; do
    echo $i
    # replace paths as required: super cluster k-mer database, sample database (plastids and mitos removed), output name 
    kmc_tools simple /disk2/hbecher_tmp/GSproject/DBsRE/DBsUnique/su$i.db21 /disk2/hbecher_tmp/GSproject/DBsFastp/$s"noPlMt" intersect $s/intersect$s"RE"$i -ocright
  done
done
