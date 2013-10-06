% Input - reads of constant length
% output:
% uniqueReads - list of unique reads packed as 64bit words
% uniqueReads_length - the number of times each read appears
% varargout is uniqueReads_inds - the original indices each unique read
%
% TEMP! added debug printings 
% 
function [uniqueReads,uniqueReads_length,packed_seqs, varargout]=createUniqueReadsAndPack(reads,readLength)

% check is reads are char array - pack them to 64bit words
if ischar(reads)
  reads_are_char = 999
  [packed_seqs, seqs_len] = pack_seqs(reads, 64);
else % already packed
  reads_are_ints = 11
  packed_seqs = reads;
end

[uniqueReads,uniqueReads_inds] = extract_sub_kmers(packed_seqs, readLength*ones(size(reads,1),1),readLength, 1,0);
clear reads
% Each read is packed as a 64bit word

% find the number of appearances of each read
[~, ~, uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));% uniqueReads_length is number of appearances of each read

varargout{1} = uniqueReads_inds;
