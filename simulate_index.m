function [i] = simulate_index(w,sum_w)
% Simulate index of the interval from which we want to sample the proposal
i=1;

u=rand(1,1);
cum_w=w(1)/sum_w;

while (u>cum_w)
   i=i+1;
   cum_w=cum_w+w(i)/sum_w;
end

end

