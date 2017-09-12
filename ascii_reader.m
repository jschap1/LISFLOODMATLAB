function [dem, ncols, nrows, xllcorner, yllcorner, cellsize] = ascii_reader (filename) 

% [DEM, ncols, nrows, xllcorner, yllcorner, cellsize] = ascii_reader (filename) 

% This function reads arc .asc files
% requires filename string as input

% j neal
% 18/6/07

%% read an ascii file
fin = fopen(filename,'r');
A = fscanf(fin,'%s',1); ncols = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); nrows = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); xllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); yllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); cellsize = fscanf(fin,'%f',1);        %#ok<NASGU>
A = fscanf(fin,'%s',1); nodata = fscanf(fin,'%f',1);          %#ok<NASGU>
dem = fscanf(fin,'%f',[ncols, nrows]);
dem = dem';

fclose('all');