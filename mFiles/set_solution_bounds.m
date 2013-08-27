% Find bounds on the value of x_i in any possible solution to an
% over-determined system Ax=y, with the added positivity constraint x_i>=0
%
% Input:
% A - coefficient matrix
% x - the correct solution
% y - free column
% sum_to_one_flag - should the solution sum to one
% run_lin_prog - should be 0 if only calculates the dependent, 1 if tries to solve
% epsilon - numerical precision (anything less is considered zero)
%
% Output:
% unique_inds - indices i's for which x_i is determined uniquely
% x_min - minimum value of each x_i in any solution
% x_max - maximum value of each x_i in any solution
%
function [unique_inds x_min x_max] = set_solution_bounds(A, x, y, sum_to_one_flag, run_lin_prog,epsilon)

if(~exist('sum_to_one_flag', 'var') || isempty(sum_to_one_flag))
    sum_to_one_flag = 0;
end
if(~exist('epsilon', 'var') || isempty(epsilon)) % set precision
    epsilon = 10^(-9);
end
disp(['uses epsilon=',num2str(epsilon),' ; change if needed. set_solution_bounds'])

[m n] = size(A);
if(sum_to_one_flag) % Just add a row of ones and run again
    [unique_inds x_min x_max] = set_solution_bounds([A' ones(n,1)]', x, [y' 1]', 0, run_lin_prog, epsilon);
    return;
end
null_basis = null(A); % get the null space
if ~isempty(null_basis)
  null_dim = size(null_basis,2);
  unique_inds = find( max(abs(null_basis),[],2) < epsilon );
else
  unique_inds = [1:n]';
end
x_min = zeros(n,1); x_max = ones(n,1); % Set trivial bounds
x_min(unique_inds) = x(unique_inds); x_max(unique_inds) = x(unique_inds);  % set unique indices


lin_prog_solver = 'cvx'; 
if(run_lin_prog)
  keyboard  
  non_unique_inds = setdiff(1:n, unique_inds); % loop on coordinates not uniquely determined
    disp(['length non_unique_inds: ',num2str(length(non_unique_inds))])
  for i=1:length(non_unique_inds)
        run_nonunique_ind = i
        f = zeros(n,1); f(non_unique_inds(i)) = 1; % set function to minimize s_i
        g =  null_basis(non_unique_inds(i),:)'; % minimize s_i. g gives the coefficients of s_i in null basis 
%keyboard
        switch lin_prog_solver 
          
            case 'matlab'
                s_min = linprog(g,-null_basis,x); % ,A,zeros(m,1));  % Solve a linear program: minimize s_i - x_i        
                s_max = linprog(-g,-null_basis,x);
            case 'cvx' % use cvx to solve 

             
             cvx_begin
               cvx_quiet(true)    
               variable s_min(null_dim)
                minimize(  g'*s_min  );
                subject to
                (-null_basis)*s_min <= x;
             cvx_end

             cvx_begin
                cvx_quiet(true)   
                variable s_max(null_dim)
                minimize(  -g'*s_max  );
                subject to
                (-null_basis)*s_max <= x;
             cvx_end
             
             
        end
       
                
        x_min(non_unique_inds(i)) = x(non_unique_inds(i)) + sum(s_min.*g);
        %    s_max = linprog(-f,-eye(n),x,A,zeros(m,1));  % Solve a linear program: maximize s_i - x_i
        x_max(non_unique_inds(i)) = x(non_unique_inds(i)) + sum(s_max.*g);
        
        if( abs(x_min(non_unique_inds(i)) -x_max(non_unique_inds(i))) < epsilon ) % add as a uniquely determined index
            
          unique_inds = [unique_inds; non_unique_inds(i)];
        end
    end
end % run linear programming to get better bounds




