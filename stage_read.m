function [metadata, stagedata] = stage_read (filename, num_locations, checkpoint)

% This function reads lisflood .stage files
% [metadata, stagedata] = stage_read (filename, num_locations)
%
% if checkpointing was used an additional argument is required to look for
%  checkpointed data
% [metadata, stagedata] = stage_read (filename, num_locations,
% checkpoint=1);
%
% The code will also check for repetition and remove duplicate lines
%
% J Neal
% 26/6/09

if nargin < 2, 
   error('Requires at least two input arguments'); 
end
if nargin < 3, 
   checkpoint = 0;
end

fin = fopen(filename, 'r');

header1 = fgets(fin);  %#ok<NASGU>
header2 = fgets(fin);  %#ok<NASGU>
header3 = fgets(fin);  %#ok<NASGU>
% reads the metadata
metadata = fscanf (fin, '%f', [4, num_locations]);
metadata = metadata';

header4 = fgets(fin);  %#ok<NASGU>
header5 = fgets(fin);  %#ok<NASGU>
header6 = fgets(fin);  %#ok<NASGU>
header7 = fgets(fin);  %#ok<NASGU>

if checkpoint == 0
% reads stage data
stagedata = fscanf (fin, '%f', [num_locations+1,inf]);
stagedata = stagedata';
else
  % scan for checkpointing lines
  line = 1;
  checks = 0;
  checkpoint = 0;
  while 1
      tline = fgetl(fin);
      % this will return -1 at end of file
      if tline == -1
          checks = checks+1;
          checkpoint(checks+1) = line; %#ok<AGROW>
          line = line-1;
          break;
      end
      TF = strcmp (tline, '####################################################### Checkpoint restart ########################################################');
      if TF == 1
          % this line is a checkpoint line
          checks = checks+1;
          checkpoint(checks+1) = line; %#ok<AGROW>
      end
      % the data could be read in at this point but this could be very
      % slow for large stage files
      line = line+1;
  end
  %now we know how may times it has checkpointed and where
  disp (['Number of checkpoints = ',num2str(checks)]);
  fclose('all');
   
  % re open file and read header again!
  fin = fopen(filename, 'r');
  header1 = fgets(fin);  %#ok<NASGU>
  header2 = fgets(fin);  %#ok<NASGU>
  header3 = fgets(fin);  %#ok<NASGU>
  % reads the metadata
  metadata = fscanf (fin, '%f', [4, num_locations]);
  metadata = metadata';
  header4 = fgets(fin);  %#ok<NASGU>
  header5 = fgets(fin);  %#ok<NASGU>
  header6 = fgets(fin);  %#ok<NASGU>
  header7 = fgets(fin);  %#ok<NASGU>
   
  % read stage data
  if checks == 2
      % reads stage data as there are no checkpoints
      stagedata = fscanf (fin, '%f', [num_locations+1,line]);
      stagedata = stagedata';
  else
      stagedata = [];
      for i = 2:checks
          tempstagedata = fscanf (fin, '%f', [num_locations+1,checkpoint(i-1)-checkpoint(i)-1]);
          tempstagedata = tempstagedata';
          % this is not efficient but shjould only happen a few time
          stagedata = [stagdata;tempstagedata];
          % read the checkpoint line
          tline = fgetl(fin); %#ok<NASGU>
      end
      % we now have lots of data that may overlap so we only need unique
      % rows
      stagedata = unique(stagedata,'rows');
  end
       
end
fclose('all');