% calculates weighted sensitivity and weighted specificity as in the manuscript
% input:
% correct_mixture_indices: indices of true mixture bacteria
% correct_mixture_freq: frequency of true mixture bacteria
% COMPASS_indices: indices of COMPASS bacteria 
% COMPASS_freq: frequency of COMPASS bacteria
% userDir,setParameters, as in other functions

% output:
% SequenceSimilarityThresholds: a vector of sequence similarity values - from 0 to 5%, namely MM0% to MM5%
% weightedSensitivity: weighted sensitivity for different sequence similarity thresholds
% weightedSpecificity: weighted specificity for different sequence similarity thresholds

% for an example see: example_evaluating_simulation_results.m

function [weightedSensitivity,weightedSpecificity,SequenceSimilarityThresholds]=calc_weighted_sensitivity_specificity(correct_mixture_indices,correct_mixture_freq,COMPASS_indices,COMPASS_freq,userDir,setParameters)


% create a fasta file of sequences that appear either in the correct or were found in by COMPASS 
if size(COMPASS_indices,1)>size(COMPASS_indices,2)
  COMPASS_indices = COMPASS_indices';
end
tmpInd = unique([correct_mixture_indices,COMPASS_indices]);

%%%%%%%%%%%%%%%%%%5
% write temporary fasta file of these sequences to be applied by mothur

% load the Sequences for the current block
load(setParameters.basicSeqKey,'positionInPart','len_uni')
numBAC = length(tmpInd); % number of sequences per block

Sequence1 = cell(numBAC,1);
for i=1:numBAC
  clear seq_*
  load([setParameters.basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
    
end
clear seq_* positionInPart
Sequence1_char = cell(numBAC,1);
Header = cell(numBAC,1);
for i=1:size(Sequence1,1)
  Sequence1_char{i} = int2nt(unpack_seqs(Sequence1{i},len_uni(tmpInd(i)),64));
  Header{i} = num2str(i); % needed for fastawrite
end

fileName = ['tmpFastaToBeDeleted_',num2str(num2str(randi(1E8)))];
fastawrite([userDir,'results/tmpDir/',fileName,'.fa'],Header,Sequence1_char)
% end write temporary fasta file
%%%%%%%%%%%%%%%%%%5



%%%%%%%%%%%%%%%%%%%%
% run mothur to calculate distances 
disp('starting mothur');
%keyboard

PWD = pwd;
cd([userDir,'results/tmpDir']);

if isunix
  w = ['!',userDir,'mothur/Linux64/mothur/mothur "#align.seqs(candidate=',...
     fileName,'.fa'...
     ', template=',userDir,'mothur/core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
     fileName,'.align, output=square)"'];
elseif ispc
  w = ['!',userDir,'mothur\Win64\mothur\mothur "#align.seqs(candidate=',...
     fileName,'.fa'...
     ', template=',userDir,'mothur\core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
     fileName,'.align, output=square)"'];
  
elseif ismac
  w = ['!',userDir,'mothur/MAC64/mothur/mothur "#align.seqs(candidate=',...
     fileName,'.fa'...
     ', template=',userDir,'mothur/core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
     fileName,'.align, output=square)"'];
  
end

eval(w)

disp('Finished calculating distances');

% distmat - the distance matrix of all the sequences
distmat=dlmread([fileName,'.square.dist'],'\t',0,1);
distmat=distmat(2:end,1:end-1);

delete([fileName,'.align.report'])
delete([fileName,'.align'])
delete([fileName,'.square.dist'])
delete(['mothur.*.logfile'])
delete([fileName,'.fa'])

cd(PWD)
% end run mothur to calculate distances 
%%%%%%%%%%%%%%%%%%


% calculate weighted sensitivity and specificity (referred to as recall and precision)

% create two vectors that hold the correct and reconstructed frequencies 
correct_mixture_vec = zeros(length(tmpInd),1);
[~,i1,i2] = intersect(correct_mixture_indices,tmpInd);
correct_mixture_vec(i2) = correct_mixture_freq(i1);

COMPASS_vec = zeros(length(tmpInd),1);
[~,i1,i2] = intersect(COMPASS_indices,tmpInd);
COMPASS_vec(i2) = COMPASS_freq(i1);
frequencySimilarityThreshold_additive=0.002; % delta

[weightedSensitivity,weightedSpecificity,SequenceSimilarityThresholds]=compareTwoSolutions(distmat,correct_mixture_vec,COMPASS_vec,setParameters);
% calculate sensitivity and specifityfy






%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arec,aprec,thresh]=compareTwoSolutions(distmat,origFreq,recFreq,setParameters)
% calculate weighted recall and precision for several values of sequence similarity threshold "sequenceSimilarityThreshold"
% the frequency similarity threshold is given by "frequencySimilarityThreshold_mult" and "frequencySimilarityThreshold_additive" for the multiplicative and additive differences (see manuscript)

arec=[];
aprec=[];
thresh=[];
for sequenceSimilarityThreshold=0:0.001:0.05
    [wrec,wprec,~]=calcWeightedRecPrec(distmat,origFreq,recFreq,sequenceSimilarityThreshold,setParameters.frequencySimilarityThreshold_mult,setParameters.frequencySimilarityThreshold_additive);
    arec=[arec wrec];
    aprec=[aprec wprec];
    thresh=[thresh sequenceSimilarityThreshold];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%
function [rec,prec,freq]=calcWeightedRecPrec(distmat,origFreq,recFreq,sequenceSimilarityThreshold,frequencySimilarityThreshold_mult,frequencySimilarityThreshold_additive)
% Calculate the weighted recall/precision of the reconstruction
% input:
% origFreq - frequencies of original seqs (0 if not present)
% recFreq - frequencies of reconstructed seqs (0 if not present)
% sequenceSimilarityThreshold - the sequence similarity threshold to define sequences as identical for the reconstruction
% output:
% rec,prec - recall (fraction of original seqs) and precision (fraction of correct seqs) of the reconstruction.
% freq - the L2 distance of the frequency vectors

origFreq(origFreq==0)=-1;
recFreq(recFreq==0)=-1;
recPresent=find(recFreq>=0);
origPresent=find(origFreq>=0);

totrec=zeros(length(origPresent),1);
totfreq=zeros(length(origPresent),1);
isrecok=zeros(length(recPresent),1);

% variables for frequency threshold
isrecokfreq=zeros(length(recPresent),1);
isorigokfreq=zeros(length(origPresent),1);

%disp(['multiplicative frequency similarity threshold is ',num2str(frequencySimilarityThreshold_mult)])
%disp(['additive frequency similarity threshold is ',num2str(frequencySimilarityThreshold_mult)])


freq=0;
summis=0;
sumfound=0;
sumnokfreq=0;
numok=0;

[neardist,nearpos]=min(distmat(origPresent,recPresent));
% go over all original sequences
for cseq=1:length(origPresent)
    corig=origPresent(cseq);
    % find corresponding reconstructed seqs
    corseqs=find(nearpos==cseq);
    % select only ones similar enough
    corseqs=corseqs(neardist(corseqs)<=sequenceSimilarityThreshold);
    if (~isempty(corseqs))
        isrecok(corseqs)=1;
        ofreq=origFreq(corig);
        rfreq=sum(recFreq(recPresent(corseqs)));
        if (((ofreq/rfreq<frequencySimilarityThreshold_mult)&&(ofreq/rfreq>1/frequencySimilarityThreshold_mult))||(abs(ofreq-rfreq)<frequencySimilarityThreshold_additive))
            freq=freq+(ofreq-rfreq)^2;
            isrecokfreq(corseqs)=1;
            isorigokfreq(cseq)=1;
            sumfound=sumfound+ofreq;
            numok=numok+1;
        else
            sumnokfreq=sumnokfreq+ofreq;
            freq=freq+ofreq^2;
        end
    else
        ofreq=origFreq(corig);
        summis=summis+ofreq;
    end
    totrec(cseq)=length(corseqs);
    totfreq(cseq)=sum(recFreq(recPresent(corseqs)));
end

tprec=sum(recFreq(recPresent(isrecokfreq>0)))/sum(recFreq(recPresent));
trec=sum(origFreq(origPresent(isorigokfreq>0)));
freq=freq+sum(recFreq(recPresent(isrecokfreq==0)).^2);
rec=trec;
prec=tprec;
