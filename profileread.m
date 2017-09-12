function [data, M] = profileread (resroot, numfiles, t, string1, string2)

% This function reads LISFLOOD-FP .profile files
%
% DATA = PROFILEREAD (RESROOT)
%   where RESROOT is a string containing the resroot from the LISFLOOD-FP
%   .par file. Absolute path names can also be used e.g. 'C:\LISFLOOD\resroot'
%   DATA is a structure containing .header {cell aray}
%                                  .segmentN (numerical array for segment N)
%
% DATA = PROFILEREAD (RESROOT, NUMFILES)
%   NUMFILES = number of files to read. Reads m files assuming that the files 
%   are numbered such that resroot-000m.profile. 
%   DATA(m) is a structure containing m files
% 
% DATA = PROFILEREAD (RESROOT, NUMFILES, 1)
%   specifies the profiles were exported a user specified time such that
%   resroot-000x-t.profile
%
% [DATA, M] = PROFILEREAD (RESROOT, NUMFILES, 1, string1, string2)
%   instructs the function to plot string1 or sting1 against string2... these
%   strings muct match the headernames in the .profile file e.g. 'Flow'
%   returns movieframes in the structure M. Also plots the data.
%
% Note: LISFLOOD-FP will name multiple files 0 to m-1 but DATA(m) will count from 1 to m
%
% j.neal
% 12/5/2008
%%
if nargin < 1, 
   error('Requires filename'); 
end
if nargin < 2,
   % only one file
   numfiles = 0;
end
if nargin < 2,
   t = 0;
end
if nargin < 4,
   % no plotting required
   plotting = 0;
   M = 0;
elseif nargin < 5,
   % use plotter
   plotting = 1;
elseif nargin < 6,
   plotting = 2;
end
%% Arrange data read
if t == 1
   fileex = '-t.profile';
else
   fileex = '.profile';
end
if numfiles == 0
   % read in a single .profile file
   filename = [resroot,fileex];
   data = read_profile (filename);
elseif numfiles < 9999
   % read a series of .profile files
   for i = 1:numfiles
       ii = i-1;
       if ii < 10
           resroot2 = [resroot,'-000',num2str(ii),fileex];
       elseif ii < 100
           resroot2 = [resroot,'-00',num2str(ii),fileex];
       elseif ii < 1000
           resroot2 = [resroot,'-0',num2str(ii),fileex];
       else
           resroot2 = [resroot,'-',num2str(ii),fileex];
       end
       %returnes structure with format data(i).segmentX... also includes .header 
       [data(i), seg_num] = read_profile (resroot2); %#ok<AGROW>
   end
else
   error ('Too many files or incorrect input to numfiles'); 
end
%% Plotting under development
% account for 0 seg num
seg_num = seg_num +1;
% if plotting == 0 % do nothing... no plots
if plotting == 1
   % find the column you want to plot
   I = strmatch(string1 , data(1).header);
   % plot a single variable
   A = size(data);
   for i = 1:A(2)
       figure(i);
       for j = 1:seg_num
           k = j-1;
           subplot(1,seg_num,j);
           eval(['plot(data(i).segment',num2str(k),'(:,',num2str(I),'));']);
           tstring = ['Plot of ',resroot,num2str(i),' segment ',num2str(k)];  
           title(tstring);
           % get movie frame
           M(i,j) = getframe(gcf); %#ok<AGROW>
       end
   end
elseif plotting == 2
   % find the columns you want to plot
   I = strmatch(string1 , data(1).header);
   I2 = strmatch(string2 , data(1).header);
   A = size(data);
   for i = 1:A(2)
       figure(i);
       for j = 1:seg_num
           k = j-1;
           subplot(1,seg_num,j);
           eval(['plot(data(i).segment',num2str(k),'(:,',num2str(I),'),data(i).segment',num2str(k),'(:,',num2str(I2),'));']);
           % get movie frame
           tstring = ['Plot of ',resroot,num2str(i),' segment ',num2str(k)];  
           title(tstring);
           M(i,j) = getframe(gcf); %#ok<AGROW>
       end
   end
end
%% read .profile file
function [data, seg_num] = read_profile (filename) 
fin = fopen(filename,'r');
more_segs = 1;
% read 'Channel_segment:' from file
temp = fscanf(fin,'%s',1); %#ok<NASGU>
while more_segs == 1
   % read filename
   seg_num = fscanf(fin,'%f',1);
   for i = 1:11
       data.header{i} = fscanf(fin,'%s',1);
   end
   testline = 1;
   j = 0;
   while testline == 1
       j = j+1;
       % reads the first bit of the line as a string
       A = fscanf(fin, '%s',1);
       % this string is either a numer of a new channel sgment of the end
       % of the file
       % if end of file
       ST = feof (fin);
       if ST == 1
           % reached the end of the file exit both while loops
           testline = 0;
           more_segs = 0;
       else
           TF = strcmp('Channel_segment:', A);
           if TF == 1
               % new segment has been found
               testline = 0;
           else
               % read rest for line into data
               B(j,1) = str2double(A); %#ok<AGROW>
               for k = 2:11
                   B(j,k) = fscanf(fin, '%f',1); %#ok<AGROW>
                end
            end
        end
    end
    eval(['data.segment',num2str(seg_num),' = B;']);
    clear B
end
fclose('all');