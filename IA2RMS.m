function [x]=IA2RMS(f,S,M,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent Adaptive^2 Rejection Metropolis Sampling (IA2RMS)
% f = target density (known up to the normalisation constant)
% for example, f = @(x) exp( -x.^2/2).*(1+(sin(3*x)).^2).*(1+(cos(5*x).^2));
% S = initial support points (at least 2, use more to avoid numerical problems) 
% M = number of samples we required
% type= 0/1/2 construction proposal (0: uniform intervals; 1: linear
% intervals; 2: exponential intervals (linear in log-scale);)
% (for heavy tailed target pdf, it is better to change the construction of the tails)
% Example with one tail:
%%>> f = @(x) exp(-x).* (x > 0); 
%%>> S=[-1 0 1]; 
%%>> x=IA2RMS(f,S,1000,0); 
% For showing results (autocorrelation etc.): showRes=1;
% For showing some interesting plot: showPlot=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results and plot (1) or not (0)
showRes=0;
showPlot=0;

% Make sure at least one point s' gives the target f(s')>0
% ( just to avoid numerical problems )
% we need at least one support point that is not in the tails
f=@(x)10^7*f(x);
S = checkInitPoints(f,S);
if isempty(S)==1
    disp('There is a problem with your initial support points, please choose others.')
    x = [];
    return
end

% Always sort support points
S = sort(S);
N = length(S);

% Evaluate target at all support points
f_S=f(S);

% Simulated variables
x=zeros(1,M);

% Count number of iterations
count=0;

% Count number of simulated points;
k=1;

% BUILD PROPOSAL from a set of support points S
[m,b,area,sum_area] = build_proposal(f_S,S,N,k,type);

% SAMPLING AND EVALUATE PROPOSAL: draw x(1) from and evaluate the proposal
[x_prev,fp_prev,j_prev]=Sample_Eval_Proposal(m,b,area,sum_area,S,N,type);
x(1)=x_prev;
ft_prev=f(x_prev);

while k<=M 

	% SAMPLING AND EVALUATE POPOSAL: draw x'from and evaluate the proposal
    [x_prop,fp_now,j_now]=Sample_Eval_Proposal(m,b,area,sum_area,S,N,type);
    
    % EVALUATE x'
    ft_now=f(x_prop);
    
    alpha1=ft_now/fp_now;
    
	if rand() > alpha1  %%%% RS test
        [f_S,S,m,b,area,sum_area] = update_proposal(x_prop,ft_now,f_S,S,m,b,area,j_now,N,k,type);
        
        % One point has been added to the set S, so N is increased
        N=N+1;
        
        % If x and x_prev are from the same interval, p_prev and j_prev have to be updated.
        if (j_prev==j_now)
            if (x_prev<x_prop)
                fp_prev=Eval_Proposal(x_prev,m(j_prev),b(j_prev),j_prev,type,N);
            else
                fp_prev=Eval_Proposal(x_prev,m(j_prev+1),b(j_prev+1),j_prev+1,type,N);
            end
        end
        
        % Since we have added one interval, if x_prop<x_prev index, j_prev has to be updated
        if (x_prop<x_prev)
            j_prev=j_prev+1;
        end
            
	else
		% accept
		q_prev=min(ft_prev,fp_prev);
		q_now=min(ft_now,fp_now);
		rho=(ft_now*q_prev)/(ft_prev*q_now);
		alpha2=min(1,rho);
        
		if rand()<=alpha2 %%%% MH test
            % accept
			x(k+1) = x_prop;
            y_aux=x(k);
            f_y=ft_prev;
            j_y=j_prev;
            
            %%%%%%%%
            %%%%%%%%
            alpha3=fp_prev/ft_prev;
            %%%%%%%%
            %%%%%%%%
            x_prev=x_prop;
            fp_prev=fp_now;
            ft_prev=ft_now;
            j_prev=j_now;
            
        else
            % reject MH test
			x(k+1) = x(k);
            y_aux=x_prop;
            f_y=ft_now;
            j_y=j_now;
            %%%%%%%%
            %%%%%%%%
            alpha3=1/alpha1;
            %%%%%%%%
            %%%%%%%%          
        end
   
    %%%% second control of IA2RMS
	if (rand()>alpha3) && (k~=M) %%%% second condition just for plotting
        [f_S,S,m,b,area,sum_area] = update_proposal(y_aux,f_y,f_S,S,m,b,area,j_y,N,k,type);
        
        % One point has been added to the set S
        N=N+1;
        
        % If y_aus and x_prev are from the same interval, p_prev and j_prev have to be updated.
        if (j_prev==j_y)
            if (x_prev<y_aux)
                fp_prev=Eval_Proposal(x_prev,m(j_prev),b(j_prev),j_prev,type,N);
            elseif (x_prev>y_aux)
                fp_prev=Eval_Proposal(x_prev,m(j_prev+1),b(j_prev+1),j_prev+1,type,N);
            end
        end
        
        % Since we have added one interval, if y_aux<x_prev index, j_prev has to be updated
        if (y_aux<x_prev)
            j_prev=j_prev+1;
        end
            
    end
		k=k+1;
    end
    count=count+1;
end
%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% RESULTS
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
if showRes==1
	% Output stats
	fprintf ('total number of iterations: %i\n', count)
	fprintf ('number of final support points: %i\n', N)
	fprintf ('linear correlation x(t),x(t+1): %f\n', sum((x(1:end-1)-mean(x)).*(x(2:end)-mean(x)))/sum((x-mean(x)).^2))
    fprintf ('linear correlation x(t),x(t+10): %f\n', sum((x(1:end-10)-mean(x)).*(x(11:end)-mean(x)))/sum((x-mean(x)).^2))
    fprintf ('linear correlation x(t),x(t+50): %f\n', sum((x(1:end-50)-mean(x)).*(x(51:end)-mean(x)))/sum((x-mean(x)).^2))
end

if showPlot==1
	% Show plots
	figure

	[bb,c]=hist(x,100);

	D=sum(bb.*(c(2)-c(1)));
	bar(c,bb*(1/D))
	hold on
    
    x1=min(x)-7:0.01:max(x)+7;
    
	fx1 = f(x1);
	fx2 = zeros(length(x1),1);
	for i=1:length(x1)
         j=sum(S<=x1(i))+1;
         fx2(i)=Eval_Proposal(x1(i),m(j),b(j),j,type,N);
	end
	D2=sum(0.01*fx1);

	plot(x1,(1/D2)*fx1,'g','LineWidth',2)
	D3=sum(0.01*fx2);
	plot(x1,(1/D3)*fx2,'r--','LineWidth',2)
    
    legend('histogram of simulated values','target distribution','proposal distribution')
end