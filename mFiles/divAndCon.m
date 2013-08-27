% function divAndCon.m: a divide and conquer wrapper
% The output if a vector X that holds the solution for each bacteria - namely it collects the results from each block
% sumRelevantReads - the number of reads thart were mapped to the matrix A - may be useful to asses the need for correcting for read errors.

function [X,sumRelevantReads]=divAndCon(uniqueReads,uniqueReads_length,curr_kp,setParameters)


readLength = setParameters.readLength;
currNumProcessors = setParameters.numProcessors;

n = length(curr_kp);

% set the group for divide and conquer: part holds the parts and tmpInd is a cell that holds the indices of bacteria in each block
if setParameters.keepOriginalOrderFlag==0
  currInd = randperm(n);
else
  currInd = 1:n;
end

if setParameters.forcePartFromOutsideFlag==0
  part = 1:setParameters.groupSize:n;
  part(end) = n+1;
else
  part = setParameters.partFromOutside;
end

tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
  tmpInd{i} = curr_kp(currInd(tmpInd{i}));
end

if length(tmpInd{end})>1.5*setParameters.groupSize
  sParts_1 = tmpInd{end}(1:setParameters.groupSize);
  sParts_2 = tmpInd{end}(setParameters.groupSize+1:length(tmpInd{end}));
  tmpInd{end} = sParts_1;
  tmpInd{end+1} = sParts_2;

  new_part_end = [part(end-1) part(end-1)+setParameters.groupSize part(end)];
  part(end-1:end) = [];
  part = [part,new_part_end];
end

% end set the group for divide and conquer
%%%%%%%%%%%%%%%5


% solve each block separately - either using the parallel toolbox or one by one
% the cell x holds the solution for each block
firstFlag = 1;
if currNumProcessors>1 
  
  % use parallel toolbox
  w = ['warning off ;matlabpool close force local;matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),[xx,ss]=doOneGroup(i,tmpInd,setParameters,uniqueReads,uniqueReads_length);x{i} = xx{i};sumRelevantReads(i) = ss(i);;end;matlabpool close;warning(''on'')'];
  eval(w);
  
else % single CPU
  
  for i=1:length(tmpInd)
    [res.x,res.sumRelevantReads]=doOneGroup(i,tmpInd,setParameters,uniqueReads,uniqueReads_length);
    x{i} = res.x{i}; % solution for block i
    sumRelevantReads(i) = res.sumRelevantReads(i); % total number of reads that were mapped to the matrix A
  end
end 
% end solving each block

if size(x{1},2)==1 
  X = zeros(1,n);
  for i=1:length(part)-1
    tmpInd = part(i):part(i+1)-1;
    tmpInd = currInd(tmpInd);
    if ~isempty(x{i})
      X(tmpInd) = x{i};
    end
  end
else
  disp('problem. divAndCon.m')
  keyboard
end



