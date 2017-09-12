function lisflood_flow_plotter (qxfile, qyfile, demfile, north, east, south, west, scaling)

% This function plots flow vectors from LISFLOOD_FP Qx and Qy files
%
% lisflood_flow_plotter (qxfile, qyfile, demfile);
%
% these input variables are strings containing the relative of full
% pathnames of the .Qx .Qy and ascii dem files to be plotted.
%
% lisflood_flow_plotter (qxfile, qyfile, demfile, N, E, S, W);
%
% N, E, W and S are optional they specify the northernmost, southeernmost
% ect cell to be displayed. this can be used to crop the domain.
%
% lisflood_flow_plotter (qxfile, qyfile, demfile, N, E, S, W, scaling);
% 
% scaling is an optional term used to manually scale the vector lengthes
% drawn by quiver (the arrow plotting function).
%
% J Neal
% 20/02/2008 

%% decide if crop is required and give imput error message
if nargin < 3, 
    error('Requires at least three input arguments'); 
end
if nargin < 7, 
    crop = 0;
else
    crop = 1;
end
if nargin < 7,
    scaling = 1;
end 
% Note don't do this for large areas ?
%% Task 1: read ascii files
% read Qx file
[QX, ncolsqx, nrowsqx, xllcornerqx, yllcornerqx, cellsizeqx] = read_file (qxfile); %#ok<NASGU>
% read Qy file
[QY, ncolsqy, nrowsqy, xllcornerqy, yllcornerqy, cellsizeqy] = read_file (qyfile); %#ok<NASGU>
% read DEM
[DEM, ncolsdem, nrowsdem, xllcornerdem, yllcornerdem, cellsizedem] = read_file (demfile);
%% crop data if required
if crop == 1;
    DEM = DEM(north:south,west:east);
    QX = QX(north:south,west:east+1);
    QY = QY(north:south+1,west:east);
    xllcornerdem = xllcornerdem + west*cellsizedem - cellsizedem;
    yllcornerdem = yllcornerdem + (nrowsdem - south) * cellsizedem;
    [nrowsqx, ncolsqx] = size(QX); [nrowsqy, ncolsqy] = size(QY); [nrowsdem, ncolsdem] = size(DEM);
end
%% Taks 2: Plot dem and generate figure
figure1 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1]);
colormap('gray');
% Create axes
axes1 = axes('Parent',figure1,'YTickLabel',{yllcornerdem + nrowsdem*cellsizedem ,yllcornerdem},...
    'YTick',[0.5,nrowsdem+0.5],...
    'YDir','reverse',...
    'XTickLabel',{xllcornerdem,xllcornerdem+ncolsdem*cellsizedem},...
    'XTick',[0.5, ncolsdem + 0.5],...
    'Layer','top',...
    'DataAspectRatio',[1 1 1]...
    %,'Clim',[13 30]... % uncomment to set colourmap limits manually
    );

box('on');
hold('all');
image(DEM,'Parent',axes1,'CDataMapping','scaled');
axis('image');
xlabel('Easting (m)');
ylabel('Northing (m)');
%% Task 3: calculate flow vectors
% pre allocate memory for loop... it will run faster this way
xlocs = zeros(ncolsqx*nrowsqx+ncolsqy*nrowsqy,1);
ylocs = xlocs;
U = xlocs;
V = xlocs;
%work out X and Y locations
for i = 1:nrowsqx
    for j = 1:ncolsqx
        xlocs((i-1)*ncolsqx+j,1) = j - 0.5;
        ylocs((i-1)*ncolsqx+j,1) = i;
        U((i-1)*ncolsqx+j,1) = QX(i,j);
    end
end
for i = 1:nrowsqy
    for j = 1:ncolsqy
        xlocs(ncolsqx*nrowsqx + (i-1)*ncolsqy +j,1) = j;
        ylocs(ncolsqx*nrowsqx + (i-1)*ncolsqy +j,1) = i - 0.5;
        V(ncolsqx*nrowsqx+(i-1)*ncolsqy+j,1) = QY(i,j);
    end
end
% work out number of nonzeros elements
num_flows = nnz(U+V);
xlocs2 = zeros(num_flows,1);
ylocs2 = xlocs2;
U2 = xlocs2;
V2 = xlocs2;
% remover zero flow locations
k = 1;
for i = 1:ncolsqx*nrowsqx+ncolsqy*nrowsqy
    if (U(i,1) == 0) && (V(i,1) == 0)
        % do nothing
    else
        xlocs2(k,1) = xlocs(i,1);
        ylocs2(k,1) = ylocs(i,1);
        U2(k,1) = U(i,1);
        V2(k,1) = V(i,1);
        k = k+1;
    end
end
% check numbers are correct
if k == num_flows + 1
    % all is ok
else
    disp('Problem with flows');
end
%% Task 4: overlay flow vectors
disp('Plotting data... if this is taking a long time you may need fewer locations');
quiver (xlocs2,ylocs2,U2,V2, scaling);
% quiver (xloc2,yloc2,U2,V2,'Parent',figure1);
disp('Done');
%% function to read LISFLOOD ascii files using filename
function [OUT, ncols, nrows, xllcorner, yllcorner, cellsize] = read_file (filename)
% read an ascii file
fin = fopen(filename,'r');
A = fscanf(fin,'%s',1); ncols = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); nrows = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); xllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); yllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); cellsize = fscanf(fin,'%f',1);        %#ok<NASGU>
A = fscanf(fin,'%s',1); nodata = fscanf(fin,'%f',1);          %#ok<NASGU>
OUT = fscanf(fin,'%f',[ncols, nrows]);
OUT = OUT';
fclose('all');