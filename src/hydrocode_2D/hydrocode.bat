@echo off

:::: Compile the program
:: msbuild hydrocode.sln /t:build /p:configuration=Debug

:::: Run the program
set EXEcute=.\hydrocode.exe
::: 
goto !

:!
goto SHARE

:SHARE
