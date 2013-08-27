% wrapper for doing several rounds of divide and conquer over the same initial set of bacteria. This is performed when the number of bacteria falls below a certain number.

% output
% currX - is the solution for the current set of bacteria

function [currX]=repeated_divAndCon(uniqueReads,uniqueReads_length,curr_kp,setParameters)

% create the blocks
if length(curr_kp)<setParameters.repeatWhenLowerThanThisValue
  setParameters.repeatRandomGroups = 2*setParameters.repeatRandomGroups;
end

numRuns = round(length(curr_kp)/setParameters.groupSize);

batchSize = round(setParameters.numBacteriaInDatabase/setParameters.groupSize);

parts = 1:batchSize:numRuns*setParameters.repeatRandomGroups;
if parts(end)~=numRuns*setParameters.repeatRandomGroups
  parts(end+1) = numRuns*setParameters.repeatRandomGroups+1;
end


clear tmpRepeatRandomGroups;
for i=2:length(parts)
  tmpRepeatRandomGroups(i) = round((parts(i)-parts(i-1))/numRuns);
end
tmpRepeatRandomGroups(1) = [];
if tmpRepeatRandomGroups(end)==0
  tmpRepeatRandomGroups(end) = 1;
end

% end creating the blocks

% run divAndCon for each partition
value_forRepeatRandomGroups = zeros(length(curr_kp),sum(tmpRepeatRandomGroups));
ind_col = 0;
for i=1:length(tmpRepeatRandomGroups)
  value_forRepeatRandomGroups(:,ind_col+1:ind_col+tmpRepeatRandomGroups(i)) = oneRepeat_wrapper(curr_kp,uniqueReads,uniqueReads_length,setParameters,tmpRepeatRandomGroups(i));  
  ind_col = ind_col+tmpRepeatRandomGroups(i);
end

mx_majority = mean(value_forRepeatRandomGroups,2);
keepMajority = find(mx_majority>setParameters.thresholdForCollectingBAC);

currX = zeros(1,length(curr_kp));
currX(keepMajority) = mx_majority(keepMajority);



