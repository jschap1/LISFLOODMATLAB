% Plot LISFLOOD-FP results
%
% Todo: add saving function

close all, clear, clc
addpath('/Users/jschapMac/Documents/Codes/LISFLOODMATLAB');
resultsdir = '/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Results_short';
cd(resultsdir);
load('LF.mat');

%% Movie of flooding over the domain

% water depth (.wd or .wdfp) OR
% water surface elevation (.elev)

resroot = 'res';
fileex = '.wd';
vartype = 'water depth (m)';
num_snaps = 50;
snapint = 30*24*3600;
dem = 'res.dem';
demCS = [0 70];
depthCS = [0 2];
framesPS = 1;
movQ = 50;

LISFLOOD_mov(resroot,fileex,vartype,num_snaps,snapint,dem,demCS,depthCS,framesPS,movQ) 

%% .mxe .max (map)

figure

subplot(2,1,1)
LF.mxe(LF.mxe==-9999) = NaN; % change nodata values to NaN
imagesc(LF.X, LF.Y, LF.mxe), colorbar
title('max WSE (m)')
xlabel('easting')
ylabel('northing')

subplot(2,1,2)
LF.mxe(LF.max==-9999) = NaN; % change nodata values to NaN
imagesc(LF.X, LF.Y, LF.max), colorbar
title('max water depth (m)')
xlabel('easting')
ylabel('northing')

%% Discharge at gauge locations (time series)

num_virtual_gauges = size(LF.discharge,2);
for k=1:num_virtual_gauges
    figure
    plot(LF.timevector, LF.discharge(:,k))
    title(['Discharge at virtual gauge ' num2str(k)]);
    xlabel('time');
    ylabel('discharge (m^3/s)');
    set(gca,'FontSize',14)
end

%% Plot boundary inflows

% NOT ENTIRELY WORKING...

% bdyname = '/Users/jschapMac/Desktop/LISFLOOD/arkansas_setup/arkansas.bdy';
% bdy = dlmread(bdyname, '\t', 4+38*4018+100, 0);

bdyname = '/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/LF_Inputs/run2/tuolumne.bdy';
bdy = dlmread(bdyname, '\t', 14000, 0);
inflow = bdy(:,1); figure
plot(inflow), title('inflow (m^2/s)')
res=1000; figure
plot(res*inflow), title('inflow (m^3/s)')

