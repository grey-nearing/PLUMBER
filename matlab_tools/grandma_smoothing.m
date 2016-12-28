function X = grandma_smoothing(X,L)

% number of data points
N = length(X);

% enough valid data?
M = length(find(isnan(X)));
if N-M<2*L; error('not enough data'); end; 

% enough good data to try gap-filling?
if M/N >0.05; error('not enough good data %d, %d.',N,M); end;

% need at least one good data to start
if isnan(X(1)); error('first data point is missing'); end;

% fill in the gaps
I = find(isnan(X));
count = 0;
while ~isempty(I);
 X(I) = X(I-1);
 count = count+1;
 if count > length(X)*0.010; error('too mand consecutive missing values (>6)'); end;
 I = find(isnan(X));
end 



