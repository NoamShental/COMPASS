function A=weightedMatrix(tmpInd,globalProfileForward,globalProfileReverse,setParameters)
% create the matrix A for the bacteria that were found while normalizing the entries by the global profile

numBAC = length(tmpInd);

% load the sequences of the relevant bacteria and convert them from 64bit words to standard sequences
load(setParameters.basicSeqKey,'positionInPart','len_uni')
Sequence1 = cell(numBAC,1);
for i=1:numBAC
  clear seq_*
  load([setParameters.basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
    
end
clear seq_* positionInPart


Sequence1_char = cell(numBAC,1);
for i=1:size(Sequence1,1)
  Sequence1_char{i} = int2nt(unpack_seqs(Sequence1{i},len_uni(tmpInd(i)),64));
end
% end load the sequences of the relevant bacteria and convert them from 64bit words to standard sequences

% create the list of k-mers in the relevant sequences
[A_not_normalized, values] = BuildMixingMatrixFromSequences(setParameters.readLength,Sequence1,len_uni(tmpInd));
k_mers = char(zeros(size(values,1),setParameters.readLength));
for i=1:size(values,1)
  k_mers(i,:) = int2nt(unpack_seqs(values(i,:),setParameters.readLength,64));
end
clear values Sequence1

% 

% create the matrix A - For each non zero entry in A_not_normalized find the position of the k_mer in the original sequence. Then use the global
% profile to set the value of A
A = zeros(size(A_not_normalized));
[px,py] = find(A_not_normalized); 
for i=1:length(px)
  k = strfind(Sequence1_char{py(i)},k_mers(px(i),:));
  forwardLocation = k;
  reverseLocation = len_uni(tmpInd(py(i)))+2-(k+setParameters.readLength);
  A(px(i),py(i)) = globalProfileForward(forwardLocation)+globalProfileReverse(reverseLocation);
  
end




