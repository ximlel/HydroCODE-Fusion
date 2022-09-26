#!/bin/bash  --login

shopt -s expand_aliases

alias preplot='~/Softwares/tecplot360ex/bin/preplot'

path=$1
echo "The first script parameter is the folder path: ./"$1

for i in `find ./${path} -name '*' ! -name '.?*' -type d`
do
    for j in `find . -name "*.tec" -type f`
    do
	preplot $j
	rm -v $j
    done
done
