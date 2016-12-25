function check_consecutive(D)

T = length(D);

for t = 2:T
 if     D(t,3) == D(t-1,3)+0.5
 elseif D(t,2) == D(t-1,2)+1                && D(t,3) == 0 && D(t-1,3) == 23.5 
 elseif D(t,1) == D(t-1,1)+1 && D(t,2) == 1 && D(t,3) == 0 && D(t-1,3) == 23.5 
 else 
  error('non-consecutive dates: D(t-1) = (%d,%d,%f) - D(t) = (%d,%d,%f)',D(t-1,:),D(t,:));
 end
end

