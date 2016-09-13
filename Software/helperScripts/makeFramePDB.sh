#!/bin/bash
cd $1
proteinFrom=$(ls | grep "proteinFrom")
proteinTo=$(ls | grep "proteinTo")
temp=$(ls | grep ".pdb" | grep -v "_frames" | grep -v "_protein" | sort -n -t _ -k 3)
output="protein_frames.pdb"
rm -f $output

# starting protein
echo "MODEL" >> $output
echo "REMARK "$proteinFrom >> $output
cat $proteinFrom >> $output
echo "ENDMDL" >> $output

for i in $temp
  do
  echo "MODEL" >> $output
  echo "REMARK "$i >> $output
  cat $i >> $output
  echo "ENDMDL" >> $output
  done
  
# target protein
echo "MODEL" >> $output
echo "REMARK "$proteinTo >> $output
cat $proteinTo >> $output
echo "ENDMDL" >> $output