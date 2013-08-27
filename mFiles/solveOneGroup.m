function [x]=solveOneGroup(normalizedBac,fracRelevantReads,metric)

numVariables = size(normalizedBac,2);
%keyboard
cvx_begin
  cvx_quiet(true)     
  variable x(numVariables)
  minimize( norm(normalizedBac*x-fracRelevantReads,metric) );
  subject to 
  x >= 0;
cvx_end

x = x./sum(x);










