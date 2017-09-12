% LISFLOOD post-processing wrapper
%
% Loads each LISFLOOD output into a Matlab structure array and saves 
% the result as a .mat file.
%
% Units can be found in LISFLOOD user manual. 
%
% Created 9/1/2017 JRS

clear, clc
addpath('/Users/jschapMac/Documents/Codes/LISFLOODMATLAB') 

%% INPUTS

% Specify input directory
cd /Users/jschapMac/Desktop/LISFLOOD/lisflood-fp/results

% Specify save location
saveloc = '/Users/jschapMac/Desktop/LISFLOOD/lisflood-fp/results';

prefix = 'res'; % "resroot" from the LFP parameter file
init_time = datetime([2000 1 1]); % Initial simulation time

%% .dem

[LF.dem.dem, ncols, nrows, xllcorner, yllcorner, cellsize] = ascii_reader ([prefix '.dem']);
LF.dem.ncols = ncols;
LF.dem.nrows = nrows;
LF.dem.xllcorner = xllcorner;
LF.dem.yllcorner = yllcorner;
LF.dem.cellsize = cellsize;

LF.info.dem = 'orig. input DEM; elevation (m)';

LF.X = xllcorner:cellsize:(xllcorner + cellsize*nrows - 1);
LF.Y = yllcorner:cellsize:(yllcorner + cellsize*ncols - 1);

%% .mass

LF.mass = dlmread([prefix, '.mass'], '\t', 1, 0);

% Extract number of massint time steps
nsteps = length(LF.mass);

LF.info.mass.Time = 'time (s) at which the data were saved';
LF.info.mass.Tstep = 'user-specified timestep/init. timestep (s)';
LF.info.mass.MinTstep = 'min. timestep used so far in simulation (s)';
LF.info.mass.NumTsteps = 'number of timesteps since the start of simulation';
LF.info.mass.Area = 'inundated area (m^2)';
LF.info.mass.Vol = 'volume of water in the domain (m^3)';
LF.info.mass.Qin = 'inflow discharge (m^3/s)';
LF.info.mass.Hds = 'water depth at downstream exit of the model domain (m)';
LF.info.mass.Qout = 'outflow discharge at downstream exit of the model domain (m^3/s)';
LF.info.mass.Qerror = 'volume error per send (m^3/s)';
LF.info.mass.Verror = 'vol. error per massint (m^3)';
LF.info.mass.RainInfEvap = 'cum. effect of infilt., evap., and rain over the simulation (10^3 m^3)';

%% .op .opelev
%% .profile

% LF.profile.info.ChanX = 'channel segment X location';
% LF.profile.info.ChanY = 'channel segment Y location';
% LF.profile.info.Chainage = '';
% LF.profile.info.Width = '';
% LF.profile.info.Mannings = '';
% LF.profile.info.Slope = '';
% LF.profile.info.BankZ = '';
% LF.profile.info.BedElev = '';
% LF.profile.info.WaterElev = '';
% LF.profile.info.WaterDepth = '';
% LF.profile.info.Flow = '';

%% .wd .elev .wdfp

LF.wd = NaN(nrows,ncols,nsteps);
LF.elev = NaN(nrows,ncols,nsteps);
LF.wdfp = NaN(nrows,ncols,nsteps);

% Note nsteps is currently defined as the number of mass intervals,
% however, the results may be written at different times, in general.

for t=0:(nsteps-1)
    if t<10
        n = ['000' num2str(t)];
    elseif t<100
        n = ['00' num2str(t)];
    elseif t<1000
        n = ['0' num2str(t)];
    else
        n = [num2str(t)]; % n = saveint number
    end
    [LF.wd(:,:,t+1), ~, ~, ~, ~, ~] = ascii_reader ([prefix '-' n '.wd']);
    [LF.elev(:,:,t+1), ~, ~, ~, ~, ~] = ascii_reader ([prefix '-' n '.elev']);
    [LF.wdfp(:,:,t+1), ~, ~, ~, ~, ~] = ascii_reader ([prefix '-' n '.wdfp']);
end

LF.info.wd = 'water depths (m) at each saveint';
LF.info.elev = 'water surface elevation (m) at each saveint';
LF.info.wdfp = 'floodplain-only water depths (m); depth of water above bankfull depth in cells containing a subgrid channel';

%% .mxe .max

[LF.mxe, ~, ~, ~, ~, ~] = ascii_reader ([prefix '.mxe']);
[LF.max, ~, ~, ~, ~, ~] = ascii_reader ([prefix '.max']);

LF.info.mxe = 'max WSE (m) for each pixel';
LF.info.max = 'max water depth (m) for each pixel';

%% .inittm .maxtm .totaltm

[LF.inittm, ~, ~, ~, ~, ~] = ascii_reader([prefix '.inittm']);
[LF.maxtm, ~, ~, ~, ~, ~] = ascii_reader([prefix '.maxtm']);
[LF.totaltm, ~, ~, ~, ~, ~] = ascii_reader([prefix '.totaltm']);

LF.info.inittm = 'time of initial inundation (hrs)';
LF.info.maxtm = 'time of max. inundation depth (hrs)';
LF.info.totaltm = 'total inundated time (hrs)';

%% .Qx .Qy .Qcx .Qcy
%% .Vx .Vy .maxVx .maxVy .maxVc .maxVcd 
%% .maxHaz
%% .QLx .QLy
%% .stage
%% .chmask
%% .segmask
%% _SGC_bedz.asc
%% SGC_bfdepth.asc
%% SGC_width.asc
%% .discharge

discharge = dlmread([prefix '.discharge'], '\t', 4, 0);
LF.discharge = discharge(:,2:end); % units are m^3/s

LF.info.discharge = 'discharge at location(s) specified in gaugefile (m^3/s)';

% Also extract time (massint times)
time = discharge(:,1);
LF.timevector = seconds(time) + init_time;

%% Save

save(fullfile(saveloc, 'LF'), 'LF')