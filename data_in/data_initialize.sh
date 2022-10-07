#!/bin/bash

MRun="$HOME/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop"
#MRun="octave"

rm -vf data_initialize.m
path=$1
echo "The first script parameter is the folder path: ./"$1

for i in `find ./${path} -name '*' ! -name '.?*' -type d`
do
    if [ -f "$i/value_start.m" ]; then
	echo "disp(\"run $i/value_start.m\")" >> data_initialize.m
	echo "run $i/value_start.m;" >> data_initialize.m
    fi
    if [ -f "$i/octave-workspace" ]; then
	rm -vf $i/octave-workspace
    fi
done

echo "data_initialize" | $MRun
rm -vf data_initialize.m
