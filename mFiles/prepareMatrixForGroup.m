% prepareMatrixForGroup.m prepare the matrix A for a specific group, i.e., block
% A is the matrix
% values are the k-mers corresponding to each row in the matrix

function [A values]=prepareMatrixForGroup(readLength,tmpInd,basicSeqNameDir,basicSeqKey)

% load the Sequences for the current block
load(basicSeqKey,'positionInPart','len_uni')

numBAC = length(tmpInd); % number of sequences per block

Sequence1 = cell(numBAC,1);
for i=1:numBAC
  clear seq_*
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
    
end
clear seq_* positionInPart

% seq is the list of possible sequences for each Sequence1
[A values] = BuildMixingMatrixFromSequences(readLength,Sequence1,len_uni(tmpInd));

