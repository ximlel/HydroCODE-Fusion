#!/bin/bash

shopt -s expand_aliases

alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
#alias MRun='octave'

for i in `find ./ -name '*' ! -name '.?*' -type d`
do
if [ -f "$i/value_start.m" ]; then
echo "run $i/value_start.m;" >> data_initialize.m 
fi
if [ -f "$i/octave-workspace" ]; then
echo "rm $i/octave-workspace" && rm $i/octave-workspace
fi
done

# echo "data_initialize" | MRun
# echo "rm data_initialize.m" && rm data_initialize.m
