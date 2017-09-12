function LISFLOOD_mov (resroot,fileex,vartype,num_snaps,snapint,dem,demCS,depthCS,framesPS,movQ) 

% LISFLOOD_mov generates movies of LISFLOOD-FP output
%
% LISFLOOD_mov (resroot,fileex,vartype,num_snaps,snapint,dem,demCS,depthCS,framesPS,movQ)
% 
% where 
% resroot is the LISFLOOD resroot (relative to m file or absolute)(string)
% fileex is the file extension '.wd' (can be changed to plot WSE and flows)
% vartype is the name of the variable being plotted e.g. 'Depth'
% num_snaps is the number of snapshot times
% snapint is the time interval between each LISFLOOD-FP snapshot (seconds)
% dem is the name of the demfile (relative or absolute)
% demCS is the range of dem hights to be plotted e.g. demCS = [0,20];
% depthCS is the range of depth values to be plotted e.g. depthCS = [0,10];
% framesPS number of frames per second (each plot is a simgle frame)
% movQ is the movie quality between 1 and 100 (100 is best)
% 
% 8/11/2007. J Neal.
%
% EXAMPLE:
% LISFLOOD_mov ('C:\Experiment1\res22','.wd','Water depth',24,10000,...
% 'C:\Experiment1\dem10m.dem.ascii',[0,60], [0,10],1,100);
%% Movie maker
% read dem filenam is filename and 1 specifies the 0.00 should not be NaN's
[DEM, ncols, nrows, xllcorner, yllcorner, cellsize] = read_file (dem); %#ok<NASGU>
% generate empty .avi file
% A = [resroot,'.avi'];
% mov = avifile(A);
% loop through snapshots
for i = 1:num_snaps
    if i < 10
        resfile = [resroot,'-000',num2str(i),fileex];
    elseif i < 100
        resfile = [resroot,'-00',num2str(i),fileex];
    elseif i < 1000
        resfile = [resroot,'-0',num2str(i),fileex];
    else
        resfile = [resroot,'-',num2str(i),fileex];
    end
    % read in deoth file for resfile i (2 specifies the 0.00 should be NaN
    [DEPTH, ncols, nrows, xllcorner, yllcorner, cellsize] = read_file (resfile);
    % plot data as figure 1
    plotter (DEPTH, DEM, nrows, ncols, cellsize, xllcorner, yllcorner,demCS,depthCS,snapint,i,vartype);
    % get movie frame
    M(i) = getframe(gcf); %#ok<AGROW>
    hold off;
%     mov = addframe(mov,M(i));
    % close the figure
    close all
end 
% % close .avi file
% mov = close(mov);
% generate empty .avi file
disp('Generating .avi movie file');
A = [resroot,fileex,'.avi'];
movie2avi(M,A,'fps',framesPS,'quality',movQ);
disp(A);
%% Read LISFLOOD ascii files using filename
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
% convert nodata to NaN;
for ii = 1:nrows
    for jj = 1:ncols
        if OUT(ii,jj) == nodata
            OUT(ii,jj) = NaN;
        elseif OUT(ii,jj) == -99.000
            % lisflood uses -99.000 as nodat so will also remove these
            OUT(ii,jj) = NaN;
        else
            % do nothing
        end
    end
end
%% Plotter data for movie slide
function plotter (DEPTH, DEM, nrows, ncols, cellsize, xllcorner, yllcorner,demCS,depthCS,snapint,i,vartype)
% number of colors in colormap
num_colors = 64;
% DEM and variable color map
cmap = [gray(num_colors); jet(num_colors)];
% scale DEM and DEPTH data to fit with new colormap
demCS2 = demCS(2)-demCS(1);
dem_scalar = num_colors/demCS2;
DEM2 = DEM*dem_scalar;
depthCS2 = depthCS(2)-depthCS(1);
depth_scalar = num_colors/depthCS2;
DEPTH2 = (DEPTH*depth_scalar)+num_colors+1;
%Plot the DEM
figure1 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1]);
colormap(cmap);
axes1 = axes('Parent',figure1,'YTickLabel',{yllcorner+(nrows*cellsize),yllcorner+cellsize},...
    'YTick', [1,nrows],...
    'YDir','reverse',...
    'XTickLabel',{xllcorner,xllcorner+(ncols*cellsize)-cellsize},...
    'XTick',[1,ncols],...
    'Layer','top');
box('on');
hold('all');
image(DEM2,'Parent',axes1);
axis('equal');axis('tight');
% Overlay the water depth (or other variable)
H = image(DEPTH2);
% code to make depth plot layer transparent
OP = ones(nrows,ncols);
for ii = 1:nrows
    for jj = 1:ncols
        if DEPTH(ii,jj) > 0
        else
            OP(ii,jj) = 0;
        end
    end
end
% make transparent. depth layer must be H
set(H,'AlphaData',OP);
% label colourbar
color_breaks = 5;
colorbreak = num_colors/color_breaks;
SCALE = zeros(color_breaks*2,1);
SCALE(1) = 0.0001;
for j = 1:color_breaks*2
    SCALE(j+1) = colorbreak*j;
end
SCALE(color_breaks+1) = num_colors+1;
VIS_SCALE = zeros(color_breaks*2,1);
VIS_SCALE(1) = demCS(1);
for j = 1:color_breaks-1
    VIS_SCALE(j+1) = ((demCS2/color_breaks)*j)+demCS(1);
end
VIS_SCALE(color_breaks+1) = depthCS(1);
for j = 1:color_breaks
    VIS_SCALE(j+1+color_breaks) = ((depthCS2/color_breaks)*j)+depthCS(1);
end
colorbar('YTick',SCALE,'YTickLabel',VIS_SCALE);
% Label axis
xlabel('BNG easting (m)');
ylabel('BNG northing (m)');
A = [vartype,' over DEM at time ',num2str((i*snapint)),'s'];
title(A);