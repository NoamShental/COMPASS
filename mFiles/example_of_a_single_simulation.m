% example of running a specific simulation - from creating the reads to applying COMPASS - this is the first introductory example.


% setting the paths
userDir = '~/CS/mFilesBAC/COMPASS/'; % ENTER THE PATH TO THE COMPASS DIRECTORY
addpath(genpath([userDir,'cvx/']))
addpath([userDir,'mFiles'])
% end setting the paths
%%%%%%%%%% 



%%%%%%%%%%%%%%%%%%%%%%%%
% setting the specific mixture
numBacteriaInDatabase = 410849; % this is fixed - the number of sequences in the Greengenes database used. In this case we use the full 16S database

% in this example - we selected 3 bacteria with uniform frequency
bacteria_in_mix = [1 10 1000]; % three bacteria in the mixture
correctWeight = zeros(1,numBacteriaInDatabase);
correctWeight(bacteria_in_mix) = 1/length(bacteria_in_mix); 
% end defining the correct mixture
%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%
% the read errror model for Illumina reads - see Supplementary method of the manuscript
ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 

p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; 

ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;
% end defining error model
%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%
% create reads 
auxData = struct;
auxData.basicSeqNameDir = [userDir,'database/full16S/datNoNonACGT/packed64/']; % database related issues - fixed for a specific database.
auxData.basicSeqKey= [userDir,'database/full16S/datNoNonACGT/keyNoNonACGT'];% database related issues - fixed for a specific database.
auxData.correctWeight = correctWeight; % the vector "x" 
auxData.Nreads = 10^5; % number of reads
auxData.readLength = 100; % read length 
auxData.seed_in = 1;% seed for randomizing reads creation
auxData.addNoiseFlag = 1; % 1 for adding noise, 0= no noise
auxData.ErrorStruct = ErrorStruct; % error model for the reads

% the following function simulates reads and outputs the unique reads and the number times each read appears. The reads are packed as 64bit words.
[readInput.uniqueReads,readInput.uniqueReads_length] = createReads_package(auxData); 

%%%%%%%5
% end creating the input reads




% setting COMPASS general parameters and specific parameters 


setParameters = struct;

% database related parameters
setParameters.numBacteriaInDatabase = 410849;
setParameters.basicSeqNameDir = [userDir,'database/full16S/datNoNonACGT/packed64/']
setParameters.basicSeqKey= [userDir,'database/full16S/datNoNonACGT/keyNoNonACGT']

% read related parameter
setParameters.readLength = 100;

% CPU related parameter
setParameters.numProcessors = 6; % number of cores

% COMPASS specific parameters 
setParameters.thresholdForCollectingBAC = 1e-3; % COMPASS threshold for 'interesting' bacteria. 
setParameters.groupSize = 1000; % COMPASS size of block
setParameters.smallestSetOfCollected = 1000;% after collecting the bacteria from each group - continue iterating id their number is smaller than smallestSetOfCollected

setParameters.repeatRandomGroups = 10; % see the following comment
setParameters.upperLimitForRepeat = 150000; % when the number of bacteria left in the divide and conquer step is lower than 150000 (value given by setParameters.upperLimitForRepeat), then  COMPASS iterates dividive and conquer 10 (value given by setParameters.repeatRandomGroups) times.

setParameters.repeatWhenLowerThanThisValue = 20000; % when the number of bacteria left in the dividive and conquer step is lower than 20000 (value given by setParameters.repeatWhenLowerThanThisValue),  then  COMPASS iterates dividive and conquer 2*10 (twice the value given by setParameters.repeatRandomGroups) times.

% end setting parameters
%%%%%%





% run COMPASS
solution = runCOMPASS(readInput,setParameters);
% solution is a vector that hold the reconstructed frequnecy of each  bacteria in the database





% look at COMPASS results
found_indices = find(solution>10^-2); % find those larger than 1%, these should be [1,10,1000] in this example

% plot a scatter plot of the correct mixture vs the solution
figure(1);clf
plot(correctWeight,solution,'.')
xlabel('correct frequency')
ylabel('COMPASS frequency')
title('COMPASS vs. correct frequency. Each dot corresponds to a sequence')


% in case one wants to know the names of the sequences - please load the database file:
data_full16S = load([userDir,'database/full16S/bac16s_full_without_ambiguous'])

% and look at the header and sequence of the found indices
data_full16S.Header_uni(found_indices)
data_full16S.Sequence_uni(found_indices)


