#!/bin/bash

path=$1
echo "The first script parameter is the folder path: ./"$1
ACK=$2
echo "The second script parameter must be: ACK"

for i in `find ./${path} -name '*' ! -name '.?*' -type d`
do
    if [[ "${ACK}" == "ACK" ]]; then
        rm -vf $i/*.dat
        rm -vf $i/*.h5
        rm -vf $i/*.tec
        rm -vf $i/*.plt
        rm -vf $i/*.vtk
    fi
done
