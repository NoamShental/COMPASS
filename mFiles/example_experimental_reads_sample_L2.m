% this script presents COMPASS analysis of larva sample L2. 
% The input are reads after preprocessing step (1), as describemd in the manucscript, namely converting each of the original 100nt reads to 11 reads of length 90nt. 


userDir = '~/CS/mFilesBAC/COMPASS/'; % ENTER THE PATH TO THE COMPASS DIRECTORY
addpath(genpath([userDir,'cvx/']))
addpath([userDir,'mFiles'])





% preprocessing step (2) as in the manuscript (preprocessing step (1) converted each 100nt reads to 11 90nt long reads)
fileName = [userDir,'experimentalReads/Drosophila_Illumina_sample_L2_after_step_1_preprocessing'];
globalFilterWinSize = 10; % 10 from each side
globalFilterStartPos = 40; % filter's start position - exlude the first 40
[uniqueReads,uniqueReads_length,uniqueReads_length_before_normalization,globalProfileForward,globalProfileReverse]=normalizeExperimentalReads(fileName,globalFilterWinSize,globalFilterStartPos);

% save the data after preprocessing step (2) for COMPASS
save([userDir,'experimentalReads/Drosophila_Illumina_sample_L2_after_step_2_preprocessing'],'uniqueReads','uniqueReads_length')

% save the global profile and the number of times each read appears prior to normalization - to be used in the postprocessing step
save([userDir,'experimentalReads/Drosophila_Illumina_sample_L2_globalProfile'],'globalProfileForward','globalProfileReverse','uniqueReads_length_before_normalization')

% end preprocessing step (2)






% run COMPASS
readInput = struct;
readInput.uniqueReads = uniqueReads;
readInput.uniqueReads_length = uniqueReads_length;


setParameters = struct;
setParameters.numBacteriaInDatabase = 212040;
setParameters.basicSeqNameDir = [userDir,'database/750/datNoNonACGT/packed64/'];
setParameters.basicSeqKey= [userDir,'database/750/datNoNonACGT/keyNoNonACGT_primers750'];
setParameters.thresholdForCollectingBAC = 1e-3;
setParameters.upperLimitForRepeat = 60000;
setParameters.groupSize = 1000; 
setParameters.repeatRandomGroups = 10;
setParameters.repeatWhenLowerThanThisValue = 100000;
setParameters.smallestSetOfCollected = 1000;
setParameters.numProcessors = 6;

setParameters.readLength = 90;
 
solution_sample_L2 = runCOMPASS(readInput,setParameters);

% end run COMPASS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postprocessing - see manuscript
% We use bacteria whose frequency is higher than 0.1% to create the matrix A. Since the matrix is small, then instead of normalizing the reads, we can normalize the entries themselves. In this case, the vector y is taken from the data before normalization (the variable was saved above: y_before_normalizing)  

% locate bacteria whose frequency is higher than 0.1%
tmpInd = find(solution_sample_L2>10^-3); 

% create the matrix A for the bacteria that were found while normalizing the entries by the global profile
A = weightedMatrix(tmpInd,globalProfileForward,globalProfileReverse,setParameters);

% calculate the original vector y, namely the vector before preprocessing step (2)
[~, values] = prepareMatrixForGroup(setParameters.readLength,tmpInd,setParameters.basicSeqNameDir,setParameters.basicSeqKey);
[y_before_normalizing,~] = mapReadsToRows(uniqueReads,uniqueReads_length_before_normalization,values);

% solve again using L1 norm instead of L2 norm (do not confuse with "sample_L1"!)
maxA = max(sum(A,1));
solution_after_postprocessing = solveOneGroup(A./maxA,y_before_normalizing/maxA,1); % both A and y were normalized by maxA for numenrical reasons

indices_of_found_bacteria = tmpInd(find(solution_after_postprocessing>10^-3)); % return to original indices

% end postprocessing
%%%%%%%%




% load the header ans sequence to list the found bacteria

data750 = load([userDir,'database/750/bac16s_primers750_without_ambiguous']);

headers_of_found_bacteria = data750.Header_750(indices_of_found_bacteria); 
% the Greengens header of found bacteria. Each cell may be a cell, whose size is larger than one in case more than one bacterium has the exact same sequence.
                                                                           
                                                                        
sequences_of_found_bacteria = data750.Sequence_750(indices_of_found_bacteria); % the 750bp sequence 



