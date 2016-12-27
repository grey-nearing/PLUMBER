function Xw = window_average(X,W)

% get dimensions
[N,D] = size(X);
assert(N > W*2);

% series length
L = floor(N/W)-1;

% init storage
Xw = zeros(L,W)./0;

% loop through dimensions of input data
for w = 1:W
 edex = w-1;
 for l = 1:L
  sdex = edex+1;
  edex = sdex+W-1;
  if length(find(isnan(X(sdex:edex))))<0.1*W
   Xw(l,w) = nanmean(X(sdex:edex));
  end
 end
end

% make sure all data is filled
%try
 assert(isempty(find(isnan(Xw),1,'first')))
%catch
% keyboard
%end

%% code below is for single chains.

%% window starting indexes
%sdex = 1:W:N-W;
%
%%initialize storage
%Xw = zeros(length(sdex),D)-9999;
%
%% loop through dimensions of input data
%for d = 1:D
% for n = 1:length(sdex)
%  Xw(n,d) = mean(X(sdex(n)+(0:W),d));
% end
%end
%
%
%
%
%keyboard
%
%assert(isempty(find(Xw(:)<-9990)))


