FOR  %%s IN (*.dscalar.nii) DO (call :ProcConvert %%s )

goto :END

:ProcConvert
set filename=%1
set outfile=%filename:.dscalar.nii=.gii%

C:\workbench\bin_windows64\wb_command -cifti-convert -to-gifti-ext %filename% .\GIFTI\%outfile%

goto :eof


:END
echo  Completed all conversions