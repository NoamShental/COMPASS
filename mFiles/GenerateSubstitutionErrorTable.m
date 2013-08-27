% Input:
% k - read length
% ErrorStruct - structure containing error model parameters
%
% Output:
% p_error_mat - a kx4x4 matrix representing error probabilities. P(i,j,k) =
% Pr(s'_i = k | s_i = j) where s is the original sequence and s' the noisy
% sequence
%
function p_error_mat = GenerateSubstitutionErrorTable(k, ErrorStruct)

p_error_mat = zeros(k,4,4); % set matrix


%p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
%    0.5    -1  0.22   0.6
%    0.35  0.15    -1  0.22
%    0.15  0.55  0.56   -1]; % , 16, length(p));

switch ErrorStruct.error_model % determine error model
    case {'exponential', 'exp'}
        alpha=log(ErrorStruct.final_error./ErrorStruct.baseline_error)./(k-1);
        p=ErrorStruct.baseline_error.*exp(alpha.*((1:k)-1));
    case 'const'
        p = ErrorStruct.baseline_error;
    case 'powerlaw'
end


% Column refer to original (true) base, rows refer to observed base.
% Order is A C G T. Model values are taken from Song's Genome Research
% bayescall paper.
%p_error_mat=reshape(
for i=1:4
    for j=1:4
        p_error_mat(:,i,j) = p .* ErrorStruct.p_one_nuc_error_mat(i,j);
        if(i==j)
            p_error_mat(:,i,j) = 1 + p_error_mat(:,i,j);
        end
    end
end


