% Input:
% A - the matrix
% x - the solution
% y - the read vector
% epsilon - numerical precision (anything less is considered zero)
%
% Output:
% unique_inds - indices i's for which x_i is determined uniquely

function [unique_inds] = findDependent(A, x, y, epsilon)

if(~exist('epsilon', 'var') || isempty(epsilon)) % set precision
    epsilon = 10^(-9);
end

[m n] = size(A);

null_basis = null(A); % get the null space
if ~isempty(null_basis)
  null_dim = size(null_basis,2);
  unique_inds = find( max(abs(null_basis),[],2) < epsilon );
else
  unique_inds = [1:n]';
end





