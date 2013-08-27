% an exanmple of evaluating results as in the manuscript

% setting the paths
userDir = '~/CS/mFilesBAC/COMPASS/'; % ENTER THE PATH TO THE COMPASS DIRECORY
addpath(genpath([userDir,'cvx/']))
addpath([userDir,'mFiles'])
% end setting the paths



%%%%%%%%%%%%%%%%%%%%%%%%%%5
% set COMPASS parameters
setParameters = struct;
setParameters.basicSeqNameDir = [userDir,'database/full16S/datNoNonACGT/packed64/'];
setParameters.basicSeqKey= [userDir,'database/full16S/datNoNonACGT/keyNoNonACGT'];
setParameters.minimal_relevant_freq = 10^-3; % for COMPASS - use bacteria whose reconstructed frequency is higher than this value
setParameters.frequencySimilarityThreshold_mult = 1.2; % multiplicative threshold for similarity between correct and reconstructed FREQUENCIES. see manuscript. 
setParameters.frequencySimilarityThreshold_additive = 0.002;% additive threshold for similarity between correct and reconstructed FREQUENCIES. see manuscript 
%%%%%%%%%%%%%%%%%%%%%%%%%%5
% end set COMPASS parameters




%%%%%%%%%%%%%%%%%%%%%%%%%
% load results 
res = load([userDir,'results/numOfReads'])

i = 1; % look at results of the first case, namely numOfReads=10^6

% find indices and frequency of correct mixture bacteria 
correct_mixture_indices = find(res.correctWeight_numOfReads{i});
correct_mixture_freq = res.correctWeight_numOfReads{i}(correct_mixture_indices);

%  find indices and frequency of COMPASS found bacteria
COMPASS_indices = find(res.solution_numOfReads{i}>setParameters.minimal_relevant_freq);
COMPASS_freq = res.solution_numOfReads{i}(COMPASS_indices);

[weightedSensitivity,weightedSpecificity,SequenceSimilarityThresholds]=calc_weighted_sensitivity_specificity(correct_mixture_indices,correct_mixture_freq,COMPASS_indices,COMPASS_freq,userDir,setParameters);

figure(1);clf
subplot(2,1,1)
plot(SequenceSimilarityThresholds,weightedSpecificity,'b')
xlabel('sequence similarity threshold')
ylabel('weighted specificity')

subplot(2,1,2)
plot(SequenceSimilarityThresholds,weightedSensitivity,'b')
xlabel('sequence similarity threshold')
ylabel('weighted sensitivity')


% note that the graph might be slightly non monotonous with the sequence similarity when looking at a specific simulation. The code is written in such a way that first sequence similarity is checked and only the does it check for frequency similarity (see manuscript). For higher sequence similarity thresholds more sequences are matched and some may fail frequency testing. 