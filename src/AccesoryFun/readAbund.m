function [vOut]=readAbund(file)
% This function access the files containing the species abundances, it
% should be a single column file with one header line.
    
fid=fopen(file);
tline=fgetl(fid); %read first line (column headers)
i=0;
while ~feof(fid) %read every row until last row reached
    i=i+1;
    tline=fgetl(fid); %read next line
    trow=sscanf(tline,'%f');
    data(i)=trow;
end
vOut=data; % take the vector, transform into numeric with cell2mat
end
