#!/bin/bash

ulimit -c unlimited
shopt -s expand_aliases

### Compile the program
# make clean
make


### Run the program
CPath=$(pwd)
MRun=`ls -l`
#alias MRun='~/Softwares/MATLAB/R2018a/bin/matlab -nojvm -nodisplay -nosplash -nodesktop'
alias MRun='octave'
EXE=./hydrocode.out  #EXEcutable program
TC=Two_Component


## RMI_Latini
:<<!
$EXE $TC/RMI_Latini/RMI_one $TC/RMI_Latini/RMI_one 1_Riemann_exact RMI_vertical
$EXE $TC/RMI_Latini/RMI_one $TC/RMI_Latini/RMI_one 2_GRP RMI_vertical
$EXE $TC/RMI_Latini/RMI_one_half $TC/RMI_Latini/RMI_one_half/ 2_GRP free
!
## A3_shell
:<<!
$EXE $TC/A3_shell/A3_shell_quarter  $TC/A3_shell/A3_shell_quarter 1_Riemann_exact Shell
$EXE $TC/A3_shell/A3_shell_quarter  $TC/A3_shell/A3_shell_quarter 2_GRP Shell
$EXE $TC/A3_shell/A3_shell_whole    $TC/A3_shell/A3_shell_whole   1_Riemann_exact free
$EXE $TC/A3_shell/A3_shell_whole    $TC/A3_shell/A3_shell_whole   2_GRP free
!
## Energy_Correct_Banks
:<<!
$EXE $TC/Energy_Correct_Banks/2D_Shock_interface_1wave/line_1281 $TC/Energy_Correct_Banks/2D_Shock_interface_1wave/line_1281 1_Riemann_exact RMI_vertical
$EXE $TC/Energy_Correct_Banks/2D_Shock_interface_1wave/line_400  $TC/Energy_Correct_Banks/2D_Shock_interface_1wave/line_400  1_Riemann_exact RMI_vertical
$EXE $TC/Energy_Correct_Banks/Bubble_He $TC/Energy_Correct_Banks/Bubble_He 1_Riemann_exact Sod
$EXE $TC/Energy_Correct_Banks/Bubble_He $TC/Energy_Correct_Banks/Bubble_He 2_GRP Sod
!
## Shock_Bubble_Quirk
:<<!
$EXE $TC/Shock_Bubble_Quirk/Bubble_He   $TC/Shock_Bubble_Quirk/Bubble_He  1_Riemann_exact Sod
$EXE $TC/Shock_Bubble_Quirk/Bubble_He   $TC/Shock_Bubble_Quirk/Bubble_He  2_GRP Sod
$EXE $TC/Shock_Bubble_Quirk/Bubble_R22  $TC/Shock_Bubble_Quirk/Bubble_R22 1_Riemann_exact Sod
$EXE $TC/Shock_Bubble_Quirk/Bubble_R22  $TC/Shock_Bubble_Quirk/Bubble_R22 2_GRP Sod
!


### gprof
# gprof -b -A -p -q $EXE gmon.out > pg
# gprof -b $EXE gmon.out | gprof2dot | dot -Tpng -o pg.png

### Valgrind
# valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
$EXE 
# gprof2dot -f callgrind -s callgrind.out | dot  -Tpng -o callgrind.png

### gcov 
# make html
