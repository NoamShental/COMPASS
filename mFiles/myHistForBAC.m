function [no,xo] = myHistForBAC(y,x)


%  Ignore NaN when computing miny and maxy.
ind = ~isnan(y);
miny = min(y(ind));
maxy = max(y(ind));
%  miny, maxy are empty only if all entries in y are NaNs.  In this case,
%  max and min would return NaN, thus we set miny and maxy accordingly.
if (isempty(miny))
  miny = NaN;
  maxy = NaN;
end  

if length(x) == 1
  if miny == maxy,
    miny = miny - floor(x/2) - 0.5; 
    maxy = maxy + ceil(x/2) - 0.5;
  end
  binwidth = (maxy - miny) ./ x;
  xx = miny + binwidth*(0:x);
  xx(length(xx)) = maxy;
  x = xx(1:length(xx)-1) + binwidth/2;
else
  xx = x(:)';
  binwidth = [diff(xx) 0];
  xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
  xx(1) = min(xx(1),miny);
  xx(end) = max(xx(end),maxy);
end

% Shift bins so the interval is ( ] instead of [ ).
xx = full(real(xx)); 
y = full(real(y)); % For compatibility
bins = xx + eps(xx);
edges = [-Inf bins];
nn = histc(y,edges,1);
edges(2:end) = xx;    % remove shift

% Combine first bin with 2nd bin and last bin with next to last bin
nn(2,:) = nn(2,:)+nn(1,:);
nn(end-1,:) = nn(end-1,:)+nn(end,:);
nn = nn(2:end-1,:);
edges(2) = [];
edges(end) = Inf;

no = nn';
xo = x;
