% for each read - find whether it can be found in the matrix A, whose reads are given by the variable values
% output - the vector fracRelevantReads, of the same length as values
% sumRelevantReads is the sum of y
function [y,sumRelevantReads]=mapReadsToRows(uniqueReads,uniqueReads_length,values)

% 
okFlag = 0;
pt = size(uniqueReads,1)-1;
while okFlag==0
  try
    part = 1:pt:size(uniqueReads,1);
    part(end) = size(uniqueReads,1)+1;
    
    i1 = [];
    i2 = [];
    for i=1:length(part)-1
      %i
      tmpInd = [part(i):part(i+1)-1]';
      [junk,i11,i22] = intersect(values,uniqueReads(tmpInd,:),'rows');
      i1 = [i1;i11];
      i2 = [i2;tmpInd(i22)];
    end
    okFlag = 1;
  catch me 
    pt = round(pt/2);
    okFlag = 0;
  end
end


y = zeros(size(values,1),1);
clear values
for i=1:length(i1)
  y(i1(i)) = uniqueReads_length(i2(i));
end

sumRelevantReads = sum(y);

