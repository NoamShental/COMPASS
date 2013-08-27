% COMPASS
% readInput is a struct that holds two fields:
% readInput.uniqueReads - unique read in 64bit packed format - The simulation that prepares the read already packs them
% readInput.uniqueReads_length - the number of times each read appears

% setParameters are all COMPASS parameters.


function solution=runCOMPASS(readInput,setParameters)

% randomize the seed - change 
rand('seed',sum(100*clock));

% initializing
setParameters.keepOriginalOrderFlag = 0; 
setParameters.forcePartFromOutsideFlag = 0;
setParameters.moveDependentToNextStageFlag = 1;
setParameters.metric = 2; 

if ~isfield(setParameters,'correctMeasurementForReadErrorsFlag')
  setParameters.correctMeasurementForReadErrorsFlag = 0; % if not stated - do not use read error correction
end

if setParameters.correctMeasurementForReadErrorsFlag==1
  disp('correcting for read errors')
end


% display field names
fn = fieldnames(setParameters);
fprintf('%%%%%%%%%%%\n');
fprintf('data:\n');
for i=1:length(fn)
  if ~strcmp(fn{i},'ErrorStruct') & ~strcmp(fn{i},'readsForCorrectReads') & ~strcmp(fn{i},'readsForCorrectReadsnumProcessors') 
    if isnumeric(setParameters.(fn{i})) 
      fprintf('%s: %d\n',fn{i},setParameters.(fn{i}));
    else
      fprintf('%s: %s\n',fn{i},setParameters.(fn{i}));
    end
  end % if ErrorStruct
end    
fprintf('end data\n')
fprintf('%%%%%%%%%%%%\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end display variables


% divide and conquer
clear store* 
k = 1;
curr_kp = 1:setParameters.numBacteriaInDatabase;
while length(curr_kp)>setParameters.smallestSetOfCollected 
  
  disp(['iteration number: ',num2str(k),'; number of relevant bacteria: ',num2str(length(curr_kp))])
  
  if setParameters.repeatRandomGroups>1
    
    if length(curr_kp)>setParameters.upperLimitForRepeat
      % running one random split. divAndCon.m
      [currX,currSumRelevantReads] = divAndCon(readInput.uniqueReads,readInput.uniqueReads_length,curr_kp,setParameters);
    else
      % repeat several times
      [currX] = repeated_divAndCon(readInput.uniqueReads,readInput.uniqueReads_length,curr_kp,setParameters);
      currSumRelevantReads = 0; % irrelevant in this case 
    end
    
    
  else % regular
    disp('running one random split. divAndCon.m ')
    [currX,currSumRelevantReads] = divAndCon(readInput.uniqueReads,readInput.uniqueReads_length,curr_kp,setParameters);
  
  end
  
  next_kp = curr_kp(find(currX>setParameters.thresholdForCollectingBAC));
  store_kp{k} = curr_kp;
  store_X{k} = currX;
  store_sumRelevantReads{k} = currSumRelevantReads;
  k = k +1;
  
  curr_kp = next_kp;
  
  if length(store_kp{k-1})==length(curr_kp)
    disp('did not reduce size, so breaks.')
    break
  end
end

if length(store_kp{k-1})~=length(curr_kp) % size was reduced
  store_kp{k} = curr_kp;
  store_X{k} = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve again and output the solution
warning off
do = length(store_kp);

setParameters.numProcessors = 1;
setParameters.moveDependentToNextStageFlag = 0;

curr_kp = store_kp{do};
setParameters.groupSize = length(curr_kp)-1;
[currX,currSumRelevantReads] = divAndCon(readInput.uniqueReads,readInput.uniqueReads_length,curr_kp,setParameters);

solution = zeros(setParameters.numBacteriaInDatabase,1);
solution(curr_kp) = currX;
warning on 






