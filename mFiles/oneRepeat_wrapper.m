function value_forRepeatRandomGroups=oneRepeat_wrapper(curr_kp,uniqueReads,uniqueReads_length,setParameters,tmpRepeatRandomGroups)

disp(['running with: ',num2str(length(curr_kp)),'; repeating runs: ',num2str(tmpRepeatRandomGroups)])
% assumes curr_kp is sorted

partTemplate = 1:setParameters.groupSize:length(curr_kp);
partTemplate(end) = length(curr_kp)+1;

curr_perm = randperm(length(curr_kp));
curr_kp_forRepeatRandomGroups = curr_kp(curr_perm);
currPart = partTemplate;

%keyboard
z = zeros(length(curr_kp),tmpRepeatRandomGroups);
z((curr_perm),1) = 1:length(curr_kp);

for l=2:tmpRepeatRandomGroups
  curr_perm = randperm(length(curr_kp));
  curr_kp_forRepeatRandomGroups = [curr_kp_forRepeatRandomGroups,curr_kp(curr_perm)];
  z((curr_perm),l) = 1+(l-1)*length(curr_kp):length(curr_kp)+(l-1)*length(curr_kp);
  currPart = [currPart,max(currPart)+partTemplate(2:end)-1];
end

% run divide and conquer - using fixed partitions 
setParameters.keepOriginalOrderFlag = 1;
setParameters.forcePartFromOutsideFlag = 1;
setParameters.partFromOutside = currPart;

[cX,cSumRelevantReads] = divAndCon(uniqueReads,uniqueReads_length,curr_kp_forRepeatRandomGroups,setParameters);

setParameters.keepOriginalOrderFlag = 0;
setParameters.forcePartFromOutsideFlag = 0;
setParameters = rmfield(setParameters,'partFromOutside');

%keyboard
value_forRepeatRandomGroups = cX(z);

if tmpRepeatRandomGroups==1
  value_forRepeatRandomGroups = value_forRepeatRandomGroups';
end

