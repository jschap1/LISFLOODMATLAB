% Convert from raster coordinates to map coordinates
% Useful for displaying bndpts, inflows

[~, R] = arcgridread('tuolumne_flowdir_ka.asc');
bndpts = [24, 142; 30, 124; 37, 162; 42, 104; 44, 194; 72, 58];
bndpts = bndpts + 1; % bc python uses 0 indexing
row = bndpts(:,1);
col = bndpts(:,2);
coords = pix2map(R, row, col);

writetable(table(coords), 'bndpts.txt');

%%

inflows = dlmread('inflows.txt');
inflows = inflows + 1;
coords = pix2map(R, inflows(:,1), inflows(:,2));
writetable(table(coords), 'inflows.txt');

%% 
% Displaying station locations for VIC routing model

[fd, R] = arcgridread('/Users/jschapMac/Desktop/Arkansas/LF_Inputs/arkansas.flowdir');
[nrow, ncol] = size(fd);

% Enter pixel locations from the stnloc file
% col = [21, 19, 24, 16, 29, 9, 13] - 1; % from left
% row = [8, 8, 7, 6, 6, 3, 5] - 1; % from bottom

col = [11 85 33 59 24 43 88 36 95 67 42 17 34 58 64 38 53 98 63 64 77 77 92 97 67 69 88 92 53 75 47 25 44 45 57 59 87 20 27] - 1; % from left
row = [41 41 40 39 34 34 32 30 28 27 27 24 19 15 14 14 7 5 5 41 32 29 17 15 6 4 5 6 30 27 38 17 26 26 26 10 24 40 39] - 1; % from bottom
row = nrow - row; % flipping (convention)

[lat, lon] = pix2latlon(R, row, col);
writetable(table([lon' lat']), '/Users/jschapMac/Desktop/Arkansas/GIS/stnlocs.txt');
