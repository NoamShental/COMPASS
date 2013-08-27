% function doOneGroup.m - solve a single block

% tmpInd is the list of all blocks
% groupNumber is the current block processed

% output: solution vector x for this block
% sumRelevantReads - the number of reads mapped to the matrix A


function [x,sumRelevantReads]=doOneGroup(groupNumber,tmpInd,setParameters,uniqueReads,uniqueReads_length)

% create the matrix for a the current block (i.e. group). A is the matrix and values lists the k-mers corresponding to each row in A.
[A values] = prepareMatrixForGroup(setParameters.readLength,tmpInd{groupNumber},setParameters.basicSeqNameDir,setParameters.basicSeqKey);
  
% correcting for read errors or just map the reads to rows
if setParameters.correctMeasurementForReadErrorsFlag 

  
  % a different mex was prepared for each read length depending on the word size - this can be easily simplified
      
  if setParameters.readLength>=32 & setParameters.readLength<64
    y = correctReads_for_SL2(values,uniqueReads,uniqueReads_length');
  elseif setParameters.readLength>=64 & setParameters.readLength<128
    y = correctReads_for_SL4(values,uniqueReads,uniqueReads_length');
  elseif setParameters.readLength>=128 & setParameters.readLength<256
    y = correctReads_for_SL8(values,uniqueReads,uniqueReads_length');
  elseif setParameters.readLength>=256 & setParameters.readLength<384
    y = correctReads_for_SL12(values,uniqueReads,uniqueReads_length');
  elseif setParameters.readLength>=384 & setParameters.readLength<512
    y = correctReads_for_SL13(values,uniqueReads,uniqueReads_length');
  else 
    disp('should prepare myCorrectReadsNew for this read length.')
  end
  
  
  sumRelevantReads(groupNumber) = sum(y); 
      
else % do not correct for read errors
  % map reads to rows of A
  [y,sumRelevantReads(groupNumber)] = mapReadsToRows(uniqueReads,uniqueReads_length,values);

  if isempty(find(y))
    disp('no reads found - check correcting for read errors')
  end
  %
end
clear values


% solve one group
if sumRelevantReads(groupNumber)>0 
  x{groupNumber}=solveOneGroup(A,y,2); 
else
  x{groupNumber} = zeros(size(A,2),1);
end 

if setParameters.moveDependentToNextStageFlag 

  B = full(A'*A);
  % move dependent to next stage
  unique_inds = findDependent(B, x{groupNumber}, y);
  oth = 1:size(B,2);
  oth(unique_inds) = [];
  x{groupNumber}(oth) = 2; % the are moved to the next stage
  %%%%%%%%%%
end  

