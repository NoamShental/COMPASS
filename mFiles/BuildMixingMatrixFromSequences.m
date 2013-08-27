% Prepare the mixing matrix from a set of sequences and read length
%
% Input:
% readLength - length of reads
% packed_seqs - sequences array (should be in packed form. 2bits per nucleotide)
% seqs_lens - length in nucleotide of each sequence 
%
% Output:
% ReadBySpeciesMat - matrix with 1 if read i appears in sequence j (in sparse format)
% kmers_packed - read sequences for each row of A
%
function [ReadBySpeciesMat kmers_packed] = BuildMixingMatrixFromSequences( ...
    readLength, packed_seqs, seqs_lens, matlab_word_size)

if(~exist('matlab_word_size', 'var'))
    matlab_word_size = 64; % set default: unix 64
end
if(matlab_word_size == 32)
    machine_word_str = 'uint32';
else
    machine_word_str = 'uint64';
end
if(iscell(packed_seqs)) % convert to an array (this is OK if lengths are not TOO different)
    num_species = length(packed_seqs);
    seq_lens_in_words = length_cell(packed_seqs);
    packed_seqs_mat = zeros(num_species, max(seq_lens_in_words), machine_word_str);
    for i=1:num_species
        packed_seqs_mat(i,1:seq_lens_in_words(i)) = packed_seqs{i};
    end
    packed_seqs = packed_seqs_mat; clear packed_seqs_mat; % save space
end

unique_flag = 1; % flag saying to return only unique kmers
[kmers_packed kmers_inds] = ...
    extract_sub_kmers(packed_seqs, seqs_lens, readLength, unique_flag, 0); % Call c-function 
ReadBySpeciesMat = sparse(kmers_inds(:,1), kmers_inds(:,2), 1); % generate a sparse reads matrix 



