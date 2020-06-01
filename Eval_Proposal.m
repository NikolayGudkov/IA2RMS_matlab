function [fp] = Eval_Proposal(x, m, b, j, type, N)
% Evaluate proposal at point x located in the interval j={1,...,N+1} with slope m and intercept b
%
%%%% Input: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%   N           number of support points
%   x           point to be injected into the set of support points
%   j           index of the interval, where the new point has to be injected
%   m           slope of the proposal in the interval j
%   b           intercept of the proposal in the interval j
%   type        type of proposal construction (0/1/2)
%
%%%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	fp          value of the proposal of the simulated value

% Obtain value of the proposal depending on the type of the piece
if (j==1|j==N+1)  % Tail pieces
    fp=exp(m*x+b);
else
    if (type~=1&type~=2) % Uniform middle piece
        fp=exp(b);
        
    elseif (type==1)     % Linear middle piece
        fp=m*x+b;
        
    elseif (type==2)     % Exponential middle piece
        fp=exp(m*x+b);
    end
end
