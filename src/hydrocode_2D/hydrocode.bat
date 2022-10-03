@echo off

:::: Compile the program
:: msbuild hydrocode.sln /t:build /p:configuration=Debug


:::: Run the program
set CPath=%~dp0
set MRun=D:/Softwares/MATLAB6p5/bin/win32/matlab.exe -nojvm -nodisplay -nosplash -nodesktop -r
set EXE=.\hydrocode.exe
set OMP_NUM_THREADS=4

::: GRP_Book
goto !
 %EXE% GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   1     EUL 33=1
 %EXE% GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2_GRP EUL 33=1
 %EXE% GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 1     EUL 33=1
 %EXE% GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2_GRP EUL 33=1
:!

::: 2D Riemnnn problem
goto !
   cd ../../data_in/two-dim/RP2D/RP2D_3_Quad
   %MRun% "value_start('INPUT'); quit;"
   cd %CPath%
 %EXE% RP2D_Positive/Config3  RP2D_Positive/Config3  2_GRP EUL
 %EXE% RP2D_Positive/Config7  RP2D_Positive/Config7  2_GRP EUL
 %EXE% RP2D_Positive/Config12 RP2D_Positive/Config12 2_GRP EUL
:!
