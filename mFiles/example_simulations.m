% simulating scenarios similar to those in Figure 3 in the manuscript.
% example simulations (a single case for each setting) was already saved in the directory results


% setting the paths
userDir = '~/CS/mFilesBAC/COMPASS/'; % ENTER THE PATH TO THE COMPASS DIRECORY
addpath(genpath([userDir,'cvx/']))
addpath([userDir,'mFiles'])
% end setting the paths

%%%%%%%%%%



% set COMPASS parameters - for definitions - see example_of_a_single_simulation.m
setParameters = struct;
setParameters.numBacteriaInDatabase = 410849;
setParameters.basicSeqNameDir = [userDir,'database/full16S/datNoNonACGT/packed64/'];
setParameters.basicSeqKey= [userDir,'database/full16S/datNoNonACGT/keyNoNonACGT'];
setParameters.thresholdForCollectingBAC = 1e-3;
setParameters.upperLimitForRepeat = 150000;
setParameters.groupSize = 1000; 
setParameters.repeatRandomGroups = 10;
setParameters.repeatWhenLowerThanThisValue = 20000;
setParameters.smallestSetOfCollected = 1000;
setParameters.numProcessors = 6;
setParameters.readLength = 100;
% end setting the parameters





% simulating scenarios, and saving both results and the original mixture. Example results appear in direcory results


% number of reads 1000k/  100k/  10k/  20k/  500k/  50k/
numOfReads = [1000  100  10  20  500  50]*1000;
for i=1:length(numOfReads)
  
  clear auxDataIn
  auxDataIn = struct;
  auxDataIn.userDir = userDir;
  auxDataIn.addNoiseFlag = 1;
  auxDataIn.numBacteriaInDatabase = 410849;
  auxDataIn.numberOfBacteriaInMix = 200;
  auxDataIn.powerLaw = 1;
  auxDataIn.readLength = 100; 
  auxDataIn.Nreads = numOfReads(i);
  [readInput.uniqueReads,readInput.uniqueReads_length,correctWeight_numOfReads{i}]=createRun_for_specific_simulation(auxDataIn);
  
  solution_numOfReads{i} = runCOMPASS(readInput,setParameters);
  
  % save both the correct and reconstructed mixtures
  save([userDir,'results/numOfReads'],'solution_numOfReads','correctWeight_numOfReads','numOfReads')
end


%  scatter plot of the results
figure(1);clf
load([userDir,'results/numOfReads'])
[~,ind] = sort(numOfReads,'ascend')
for i=1:length(numOfReads)
  subplot(2,3,i)
  plot(correctWeight_numOfReads{ind(i)},solution_numOfReads{ind(i)},'.')
  xlabel('correct frequency')
  ylabel('COMPASS frequency')
  title(['number of reads: ',num2str(numOfReads(ind(i)))])
end



% number of bacteria 10/  100/  1000/  200/  400/  600/
numOfBact = [10  100  1000  200  400  600]
for i=1:length(numOfBact)
  clear auxDataIn
  auxDataIn = struct;
  auxDataIn.userDir = userDir;
  auxDataIn.addNoiseFlag = 1; % read errors are added according to the model (see manuscript)
  auxDataIn.numBacteriaInDatabase = 410849;
  auxDataIn.numberOfBacteriaInMix = numOfBact(i);
  auxDataIn.powerLaw = 0; % frequency taken from a uniform distribution
  auxDataIn.readLength = 100; 
  auxDataIn.Nreads = 10^6; % number of reads
  [readInput.uniqueReads,readInput.uniqueReads_length,correctWeight_numOfBact{i}]=createRun_for_specific_simulation(auxDataIn);
  
  solution_numOfBact{i} = runCOMPASS(readInput,setParameters);
  
  % save both the correct and reconstructed mixtures
  save([userDir,'results/numOfBact'],'solution_numOfBact','correctWeight_numOfBact','numOfBact')
end

%  scatter plot of the results
figure(2);clf
load([userDir,'results/numOfBact'])
[~,ind] = sort(numOfBact,'ascend')
for i=1:length(numOfBact)
  subplot(2,3,i)
  plot(correctWeight_numOfBact{ind(i)},solution_numOfBact{ind(i)},'.')
  xlabel('correct frequency')
  ylabel('COMPASS frequency')
  title(['numOfBact: ',num2str(numOfBact(ind(i)))])
end


% read length 100/  150/  200/  75/ 50/ 35/
readLength = [100  150  200  75 50 35];
for i=1:length(readLength)
  clear auxDataIn
  auxDataIn = struct;
  auxDataIn.userDir = userDir;
  auxDataIn.addNoiseFlag = 1;% read errors are added according to the model (see manuscript)
  auxDataIn.numBacteriaInDatabase = 410849;
  auxDataIn.numberOfBacteriaInMix = 200;
  auxDataIn.powerLaw = 1; % frequency taken from 1/x 
  auxDataIn.readLength = readLength(i); 
  auxDataIn.Nreads = 10^6; % number of reads
  [readInput.uniqueReads,readInput.uniqueReads_length,correctWeight_readLength{i}]=createRun_for_specific_simulation(auxDataIn);
  
  setParameters.readLength = readLength(i); % in this case we need to change the values in setParameters since the read length changed
  solution_readLength{i} = runCOMPASS(readInput,setParameters);

  % save both the correct and reconstructed mixtures
  save([userDir,'results/readLength'],'solution_readLength','correctWeight_readLength','readLength')
end

%  scatter plot of the results
figure(3);clf
load([userDir,'results/readLength'])
[~,ind] = sort(readLength,'ascend')
for i=1:length(readLength)
  subplot(2,3,i)
  plot(correctWeight_readLength{ind(i)},solution_readLength{ind(i)},'.')
  xlabel('correct frequency')
  ylabel('COMPASS frequency')
  title(['read length: ',num2str(readLength(ind(i)))])
end
