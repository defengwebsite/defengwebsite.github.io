%% If it is a 64-bit  system with microsoft complier, use the following compile command:
mex -O -largeArrayDims mexeig.c 'C:\Program Files\MATLAB\R2012b\extern\lib\win64\microsoft\libmwlapack.lib' 'C:\Program Files\MATLAB\R2012b\extern\lib\win64\microsoft\libmwblas.lib';

 

%% If it is a 32-bit  system with microsoft complier, use the following compile command:   
% mex -O -largeArrayDims mexeig.c 'C:\MATLAB\R2011b\extern\lib\win32\microsoft\libmwlapack.lib' 'C:\MATLAB\R2011b\extern\lib\win32\microsoft\libmwblas.lib';
 
%% If it is a 32-bit  systme with lcc complier, use the following compile
%% command

% mex -O -largeArrayDims mexeig.c 'C:\MATLAB\R2011b\extern\lib\win32\lcc\libmwlapack.lib' 'C:\MATLAB\R2011b\extern\lib\win32\lcc\libmwblas.lib';
