function [m, b, area, sum_area] = build_proposal(f_S,S,N,t,type)
% Builds a proposal from a set of support points S
%
%%%% Input: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	S	 set of N support points
%   f_S  values of the target distribution at the set of support points
%   t    iterations on the chain
%   type type of proposal construction (0/1/2)
%
%%%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%   m           slopes of the proposal
%   b           intercepts of the proposal
%   area        area of the pieces of the proposal
%   sum_area    total area under the proposal

% Parameters for the initialization of the tails, to reduce the dependence of choice of the initial points. 
beta=0.95;	% 0<beta<1 	
tau=0.01;	% tau>0		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(1,N+1);       % intercepts
m = zeros(1,N+1);       % slopes
area = zeros(1,N+1);    % areas
xa=[-Inf S +Inf];       % extended set of support points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the pieces			        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type 0 for approximations of order 0 (constant pieces)
if (type~=1) & (type~=2)
    b(2:N)=max([log(f_S(1:N-1));log(f_S(2:N))]);
    
% Type 1 for approximations of order 1 (linear pieces)
elseif type == 1   
  m(2:N)=(f_S(1:N-1)-f_S(2:N))./(S(1:N-1)-S(2:N));
  b(2:N)=f_S(1:N-1)-m(2:N).*S(1:N-1);
  
% Type 2 for approximations of order 1 in the log-scase (exponential pieces)
elseif type == 2  
  m(2:N)=(log(f_S(1:N-1))-log(f_S(2:N)))./(S(1:N-1)-S(2:N));
  b(2:N)=log(f_S(1:N-1))-m(2:N).*S(1:N-1);   
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the tail pieces 			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The tails connect with the end pieces
v(1) = log(f_S(1));
v(2) = log(f_S(2));
v(3) = log(f_S(end-1));
v(4) = log(f_S(end));
 
% Control to avoid numerical problems 
pos = find(v==-Inf);
if isempty(pos) == 0
	v(pos) = -50;
end
 
% first slope and intercept
m(1) = ( v(2) - v(1) ) / ( S(2) - S(1) );
m(1)= m(1)*(1-beta*exp(-tau*t));
b(1) = v(1) - m(1) * S(1);

% end slope and intercept
m(N+1) = ( v(end) - v(end-1) ) / ( S(end) - S(end-1));
m(N+1)= m(N+1)*(1-beta*exp(-tau*t));
b(N+1) = v(end) - m(N+1) * S(end);

% Control for 'suitable' slopes, to get a proper proposal pdf  
if m(1) <= 0
	m(1) = 0.05;
	b(1) = v(1) - m(1) * S(1);
end
if m(N+1) >= 0
	m(N+1) = -0.05;
	b(N+1) = v(4) - m(N+1) * S(end);
end

%%%%%%%%%%%% COMPUTE AREA (tail pieces) %%%%%%%%%%%%%%
   area(1)=1/m(1)*(exp(m(1)*xa(2)+b(1))-exp(m(1)*xa(1)+b(1)));
   area(N+1)=1/m(N+1)*(exp(m(N+1)*xa(N+2)+b(N+1))-exp(m(N+1)*xa(N+1)+b(N+1)));

%%%%%%%%%%%% COMPUTE AREA (middle pieces) %%%%%%%%%%%
if (type~=1) & (type~=2)
    area(2:N) = (xa(3:N+1)-xa(2:N)).*exp(b(2:N)); 
    
%%%%%%%%%%%%%%%%%
elseif type==1
%%%%%%%%%%%%%%%%%    
   fs1=m(2:N).*xa(2:N)+b(2:N);
   fs2=m(2:N).*xa(3:N+1)+b(2:N);  
   area(2:N) = (fs1+fs2).*(xa(3:N+1)-xa(2:N))/2;
   
%%%%%%%%%%%%%%%%%
elseif type==2
%%%%%%%%%%%%%%%%%
   fs1=m(2:N).*xa(2:N)+b(2:N);
   fs2=m(2:N).*xa(3:N+1)+b(2:N);  
   area(2:N) = (exp(fs2)-exp(fs1))./m(2:N);
end

sum_area=sum(area);
