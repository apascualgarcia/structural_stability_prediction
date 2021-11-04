function [ gammaOut ] =readFortranGammaOut( file,Sa,Sp )
%readFortranGammaOut Reads the output file of BeyondMeanField.f
% representing the gamma matrix (long format) and transforms it into
% a wide format matrix, it requires the number of animal and plant species
% for the loop

data=readtable(file,'FileType','text','ReadVariableNames',false);
data=table2array(data); % Convert data to double
if(size(data,1) ~= Sa*Sp)
    '** Error: Size of the file do not match the expected dimensions --- abort execution'
    quit
end
gammaOut=zeros(Sa,Sp);
i=1;
j=0;
k=0;
while(k < size(data,1)) %read every row until last row reached
    k=k+1;
    j=j+1;
    gammaOut(i,j)=data(k,1); % no need to initialize a cell
    if(j == Sp) % Create a wide format matrix with the expected dimension
        i=i+1;
        j=0;
    end
end

end

