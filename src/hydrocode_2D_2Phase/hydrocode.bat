@echo off

:::: Compile the program
:: msbuild hydrocode.sln /t:build /p:configuration=Debug


:::: Run the program
set CPath=%~dp0
set MRun=D:/Softwares/MATLAB6p5/bin/win32/matlab.exe -nojvm -nodisplay -nosplash -nodesktop -r
set EXEcute=.\hydrocode.exe
set OMP_NUM_THREADS=4

::: 2D Riemnnn problem
goto !
   cd ../../data_in/two-dim/
   %MRun% "value_start; quit;"
   cd %CPath%
 %EXEcute% 
:!
