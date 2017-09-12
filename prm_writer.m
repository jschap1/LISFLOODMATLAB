function prm_writer (roughness) 

% This function writes LISFLOOD-FP parameter files numbered 1 to n using
% the roughness values in the n-by-2 matrix roughness. You will need to
% edit the code for your application.

[num_sims] = size(roughness);
for i = 1:num_sims(1)
a = ['mymodel', num2str(i),'.par']; b = num2str(i);
fout = fopen(a,'w');
fprintf(fout,'DEMfile                  mydem.dem.ascii\n');
fprintf(fout,'resroot                  myres%s\n',b);
fprintf(fout,'dirroot                  mydir\n');
fprintf(fout,'sim_time                 245000.0\n');
fprintf(fout,'initial_tstep            10.0\n');
fprintf(fout,'massint                  300.0\n');
fprintf(fout,'saveint                  10000.0\n');
%fprintf(fout,'overpass                 135000.0\n');
fprintf(fout,'fpfric                   %1.3f\n',roughness(i,2));
fprintf(fout,'nch                      %1.3f\n',roughness(i,1));
%fprintf(fout,'manningfile              mymanfile%s.n.ascii\n',b);
fprintf(fout,'riverfile                myriver.river\n');
fprintf(fout,'bdyfile                  mybdy.bdy\n');
fprintf(fout,'stagefile                mystage.stage\n');
fprintf(fout,'bcifile                  mybci.bci\n');
%fprintf(fout,'startfile                mystart.old\n');
%fprintf(fout,'checkpoint               2\n');
%fprintf(fout,'checkfile                mycheck.chkpnt\n');
fprintf(fout,'elevoff\n');
fprintf(fout,'qoutput\n');
fprintf(fout,'diffusive\n');
fclose('all');
end