function ascii_write (filename, DEM, xllcorner, yllcorner, cellsize, nodata)

% ascii_write (filename, DEM, xllcorner, yllcorner, cellsize, nodata)
% only filename (string) and DEM (2D matrix) are essential
% this function writes ascii raster files for arc
% j neal
% 6/3/2008

if nargin < 2, 
    error('Requires at least two input arguments'); 
end
if nargin < 4, 
    xllcorner = 0;
    yllcorner = 0;
end
if nargin < 5,
    cellsize = 1;
end
if nargin < 6
    nodata = -9999;
end

A = size(DEM);
fout = fopen (filename,'w');
fprintf(fout, 'ncols         %.0f\n', A(2));  
fprintf(fout, 'nrows         %.0f\n', A(1));
fprintf(fout, 'xllcorner     %f\n', xllcorner);
fprintf(fout, 'yllcorner     %f\n', yllcorner);
fprintf(fout, 'cellsize      %f\n', cellsize);
fprintf(fout, 'NODATA_value  %f\n', nodata); 
for ii = 1:A(1)
    B = DEM(ii,:);
    fprintf(fout, '%1.3f ', B);
    fprintf(fout, '\n');
end

fclose('all');