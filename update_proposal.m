function [f_S,S,m,b,area,sum_area] = update_proposal(x, f_x, f_S_old, S_old, m_old, b_old, area_old, j, N, t, type)
% UPDATE_PROPOSAL updates the proposal distribution by inserting point x in the interval j
%
%%%% Input: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	S_old       set of N support point before x has been inserted
%   f_S_old     values of the target distribution at points S_old
%   t           iterations on the chain
%   m_old       "old" slopes of the linear pieces
%   b_old       "old" intercepts of the linear pieces
%   w_old       "old" weight of the linear pieces
%   area_old    "old" area of the pieces of the proposal
%   x           point to be injected into the set of support points
%   j           index of the interval, where the new point has to be injected
%   type        type of proposal construction (0/1/2)
%
%%%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	S           set of N+1 support point after x has been inserted
%   f_S         values of the target distribution at points S after x has been inserted
%   m           "updated" slopes of the proposal
%   b           "updated" intercepts of the proposal
%   area        "updated" area of the pieces of the proposal
%   sum_area    total area under the proposal


% Parameters for the initialization of the tails, to reduce the dependence of choice of the initial points. 
beta=0.95;	% 0<beta<1 	
tau=0.01;	% tau>0	

%%%%%%%%%%%%%%%%%%%%%%
% Initiate variables %
%%%%%%%%%%%%%%%%%%%%%%

m=zeros(1,N+2);
b=zeros(1,N+2);
area=zeros(1,N+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insert x into the support grid  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=[S_old(1:j-1),x,S_old(j:N)];
f_S=[f_S_old(1:j-1),f_x,f_S_old(j:N)];
xa=[-Inf S +Inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the pieces	j: [S_j, x) and j+1: [x S_{j+1})  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type 0 for approximations of order 0 (constant pieces)
if (type ~= 1) & (type~=2)
    
    if(j==1)
        b(2:N+1)=[max(log(f_x),log(f_S_old(1))),b_old(2:N)];
    elseif(j==N+1)
        b(2:N+1)=[b_old(2:N),max(log(f_S_old(N)),log(f_x))];
    else
        b_j=max(log(f_S_old(j-1)),log(f_x));
        b_{j+1}=max(log(f_x),log(f_S_old(j)));
        b(2:N+1)=[b_old(2:j-1),b_j,b_{j+1},b_old(j+1:N)];
    end
      
% Type 1 for approximations of order 1 (linear pieces)
elseif type == 1
  
     if(j==1)
        m(2:N+1)=[(f_S_old(1)-f_x)/(S_old(1)-x),m_old(2:N)];
        b(2:N+1)=[f_x-m(2)*x,b_old(2:N)];
    elseif(j==N+1)
        m(2:N+1)=[m_old(2:N),(f_x-f_S_old(N))/(x-S_old(N))];
        b(2:N+1)=[b_old(2:N),f_x-m(N+1)*x];
     else
        m_j=(f_x-f_S_old(j-1))/(x-S_old(j-1));
        m_{j+1}=(f_S_old(j)-f_x)/(S_old(j)-x);
        b_j=f_x-m_j*x;
        b_{j+1}=f_x-m_{j+1}*x;
        m(2:N+1)=[m_old(2:j-1),m_j,m_{j+1},m_old(j+1:N)];
        b(2:N+1)=[b_old(2:j-1),b_j,b_{j+1},b_old(j+1:N)];
     end
     
% Type 2 for approximations of order 1 in log scale (exponential pieces)
elseif type == 2
  
     if(j==1)
        m(2:N+1)=[(log(f_S_old(1))-log(f_x))/(S_old(1)-x),m_old(2:N)];
        b(2:N+1)=[log(f_x)-m(2)*x,b_old(2:N)];
    elseif(j==N+1)
        m(2:N+1)=[m_old(2:N),(log(f_x)-log(f_S_old(N)))/(x-S_old(N))];
        b(2:N+1)=[b_old(2:N),log(f_x)-m(N+1)*x];
     else
        m_j=(log(f_x)-log(f_S_old(j-1)))/(x-S_old(j-1));
        m_{j+1}=(log(f_S_old(j))-log(f_x))/(S_old(j)-x);
        b_j=log(f_x)-m_j*x;
        b_{j+1}=log(f_x)-m_{j+1}*x;
        m(2:N+1)=[m_old(2:j-1),m_j,m_{j+1},m_old(j+1:N)];
        b(2:N+1)=[b_old(2:j-1),b_j,b_{j+1},b_old(j+1:N)];
     end
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
m(N+2) = ( v(end) - v(end-1) ) / ( S(end) - S(end-1));
m(N+2)= m(N+2)*(1-beta*exp(-tau*t));
b(N+2) = v(end) - m(N+2) * S(end);

% Control for 'suitable' slopes, to get a proper proposal pdf  
if m(1) <= 0
	m(1) = 0.05;
	b(1) = v(1) - m(1) * S(1);
end
if m(N+2) >= 0
	m(N+2) = -0.05;
	b(N+2) = v(4) - m(N+2) * S(end);
end

%%%%%%%%%%%% COMPUTE AREA (tail pieces) %%%%%%%%%%%%%%
   area(1)=(exp(m(1)*xa(2)+b(1))-exp(m(1)*xa(1)+b(1)))/m(1);
   area(N+2)=(exp(m(N+2)*xa(N+3)+b(N+2))-exp(m(N+2)*xa(N+2)+b(N+2)))/m(N+2);


%%%%%%%%%%%% UPDATE AREA (middle pieces) %%%%%%%%%%%
if (type~=1) & (type~=2)
  
   if (j==1)
       area(2)=(S_old(1)-x)*exp(b(2));
       area(3:N+1)=area_old(2:N);
   elseif (j==N+1)
       area(N+1)=(x-S_old(N))*exp(b(N+1));
       area(2:N)=area_old(2:N);
   else
       area_j=(x-S_old(j-1))*exp(b_j);
       area_{j+1}=(S_old(j)-x)*exp(b_{j+1});
       area(2:N+1)=[area_old(2:j-1),area_j,area_{j+1},area_old(j+1:N)];
   end   
   
%%%%%%%%%%%%%%%%%%
elseif type==1
%%%%%%%%%%%%%%%%%%
   
   if (j==1)
       area(2)=((S_old(1)+x)*m(2)+2*b(2))*(S_old(1)-x)/2;
       area(3:N+1)=area_old(2:N);  
   elseif (j==N+1)
       area(2:N)=area_old(2:N);
       area(N+1)=((S_old(N)+x)*m(N+1)+2*b(N+1))*(x-S_old(N))/2;    
   else
       area_j=((S_old(j-1)+x)*m(j)+2*b(j))*(x-S_old(j-1))/2;
       area_{j+1}=((S_old(j)+x)*m(j+1)+2*b(j+1))*(S_old(j)-x)/2;
       area(2:N+1)=[area_old(2:j-1),area_j,area_{j+1},area_old(j+1:N)]; 
   end
   
%%%%%%%%%%%%%%%%%%
elseif type==2
%%%%%%%%%%%%%%%%%%
   
   if (j==1)
       area(2)=(exp(S_old(1)*m(2)+b(2))-exp(x*m(2)+b(2)))/m(2);
       area(3:N+1)=area_old(2:N);  
   elseif (j==N+1)
       area(2:N)=area_old(2:N);
       area(N+1)=(exp(x*m(N+1)+b(N+1))-exp(S_old(N)*m(N+1)+b(N+1)))/m(N+1);
   else
       area_j=(exp(x*m(j)+b(j))-exp(S_old(j-1)*m(j)+b(j)))/m(j);
       area_{j+1}=(exp(S_old(j)*m(j+1)+b(j+1))-exp(x*m(j+1)+b(j+1)))/m(j+1);
       area(2:N+1)=[area_old(2:j-1),area_j,area_{j+1},area_old(j+1:N)]; 
   end

end

sum_area=sum(area);
