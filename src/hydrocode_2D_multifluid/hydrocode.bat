@echo off

:::: Compile the program
:: msbuild hydrocode.sln /t:build /p:configuration=Debug


:::: Run the program
set CPath=%~dp0
set MRun=D:/Softwares/MATLAB6p5/bin/win32/matlab.exe -nojvm -nodisplay -nosplash -nodesktop -r
set EXEcute=.\hydrocode.exe

::: GRP_Book
goto !
 %EXEcute% GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2 1     EUL 33=1
 %EXEcute% GRP_Book/6_1_Sod_10_lines   GRP_Book/6_1_Sod_10_lines   2 2_GRP EUL 33=1
 %EXEcute% GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2 1     EUL 33=1
 %EXEcute% GRP_Book/6_1_Sod_10_columns GRP_Book/6_1_Sod_10_columns 2 2_GRP EUL 33=1
:!

::: 2D Riemnnn problem
goto !
 cd ../../data_in/two-dim/RP2D/RP2D_3_Quad
 %MRun% "value_start; quit;"
 cd %CPath%
 %EXEcute% RP2D_Positive/Config3  RP2D_Positive/Config3  2 2_GRP EUL
 %EXEcute% RP2D_Positive/Config7  RP2D_Positive/Config7  2 2_GRP EUL
 %EXEcute% RP2D_Positive/Config12 RP2D_Positive/Config12 2 2_GRP EUL
:!
