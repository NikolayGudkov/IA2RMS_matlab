function [sample,fp,j]=Sample_Eval_Proposal(m,b,area,sum_area,S,N,type)
% Samples from a proposal distribution and evaluates proposal at the simulated point
%
%%%% Input: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	S           set of N support points
%   m           slopes of the proposal
%   b           intercepts of the proposal
%   area        area of the pieces of the proposal
%   sum_area    total area under the proposal
%   type        type of proposal construction (0/1/2)
%
%%%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%   sample      simulated value
%   fp          value of the proposal at the simulated value
%   j           index of the interval used for simulation of the sample


   
% Choose interval from which we simulate: an integer between 1 and N+1.
   j=simulate_index(area,sum_area);
% Sample from the chosen interval
    if (j==1)  % Sample from the left tail
        sample=S(1)+log(rand())/m(1);
        fp=exp(m(1)*sample+b(1));
        
    elseif (j==N+1)  % Sample from the right tail        
        sample=S(N)+log(rand())/m(N+1);
        fp=exp(m(N+1)*sample+b(N+1));
        
    else  % Sample from the middle piece depending on the type of the piece
        
        if (type ~= 1 & type~=2)  % Sample for a uniform piece
            sample=S(j-1)+(S(j)-S(j-1))*rand();
            fp=exp(b(j));
            
        elseif (type == 1) % Sample from a linear piece
            C=m(j)*S(j-1)^2/2+b(j)*S(j-1)+area(j)*rand();
            sample=(-b(j)+sqrt(b(j)^2+2*m(j)*C))/m(j);
            fp=m(j)*sample+b(j);
            
        elseif (type == 2) % Sample from an exponential piece
            sample=(log(exp(m(j)*S(j-1)+b(j))+m(j)*area(j)*rand())-b(j))/m(j);
            fp=exp(m(j)*sample+b(j));
        end
    end
