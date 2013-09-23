% Find the set of species which can be reconstructed uniquely 
% from short-reads data (assuming no errors model) 
% 
% Input: 
% packed_seqs - cell array with species sequences
% seqs_lens - species sequences length in nucleotids
% L - read length
% 
% Output: 
% identifiable_inds_vec - set of identifiable species 
% 
function [identifiable_inds_vec A kmers_packed] = get_identifiable_species(packed_seqs, seqs_lens, L)

AssignGeneralConstants;
num_species = length(seqs_lens); 
identifiable_inds_vec = zeros(num_species, 1); 

sprintf('Building Read-by-Species Matrix ..')
time_to_build_mixing_matrix = cputime;
[A kmers_packed] = BuildMixingMatrixFromSequences(L, packed_seqs, seqs_lens); % use word size  
size(A)
time_to_build_mixing_matrix = cputime - time_to_build_mixing_matrix
% % % % time_build = cputime; 
% % % % [A_small small_kmers_packed] = BuildMixingMatrixFromSequences(L-20, packed_seqs, seqs_lens, 32);
% % % % time_build = cputime - time_build
% % % % time_project = cputime;
% % % % [A_projected projected_kmers_pakced] = ProjectMixingMatrixOnLowerReadLength(A, kmers_packed, L, L-20, packed_seqs, seqs_lens);
% % % % time_project = cputime - time_project
% % % % [intersect_kmers_packed I J] = intersect(small_kmers_packed, projected_kmers_pakced, 'rows'); 
% % % % fff = find(A_small(I,:) ~= A_projected(J,:))
% % % % 
XXX = 12341234
% return; 

run_lin_prog = 0; sum_to_one_flag = 1; 
x = zeros(num_species,1); y = []; epsilon = []; % fill space with zeros
sprintf('Finding Identifiable Species using Null-Space of A ..')

if((machine == PC) && (num_species < 1000)) % condition for PC
    A = full(A);  % use full matrix 
end
time_to_find_identifiable_species=cputime;
unique_inds = set_solution_bounds(A, x, y, sum_to_one_flag, run_lin_prog,epsilon);
time_to_find_identifiable_species = cputime - time_to_find_identifiable_species

identifiable_inds_vec(unique_inds) = 1; % set the identifiable species

