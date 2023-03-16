%%
clear
close all
clc

%%
a=1.5;
b=2^(a/2)*gamma(1+a/2);
C=((1-a)/(gamma(2-a)*cos(pi*a/2)))*(a~=1)+2/pi*(a==1);

% Fixed parameters
K=1e4;              % Number of summands in series representation 
n=1e4;              % sample size       
delta=0.5;          % t_i-t_{i-1}
t_eval=(0:n)*delta;
T=t_eval(end);      % Range of path, i.e. X(t), t in [0,T]

%% Set number of est frequencies N
error=1e-3;     % Decrease error for larger sample size n
N=floor((n*error)^(5/2));

%%
% Choose example
i=2;

% Optionally run loop for all four examples
% for i=1:4

% Generate seq. of r.v.'s
rng('default')  % set seed
G=normrnd(0,1,2,K);
Gam=cumsum(exprnd(1,1,K));

switch i
    case 1
        f=@(x) normpdf(x,0,1);
        Z=normrnd(0,1,1,K);
        f_range=5;
    case 2
        sample=linspace(-20,20,1e5);
        f=@(x) x.^2.*exp(-abs(x))/4;
        fsample=f(sample);
        Z=randpdf(fsample, sample, [1 K]);
        f_range=10;
    case 3
        sample=linspace(-30,30,10e5);
        f=@(x) 0.5./x.^2.*(abs(x)>=1);
        fsample=f(sample);
        Z=randpdf(fsample, sample, [1 K]);
        Z(isnan(Z))=0;
        f_range=15;
    case 4
        f=@(x) unifpdf(x,-1,1);
        Z=rand(1,K)*2-1;
        f_range=2;
end
mR=1;       % since all f are pdfs
f_eval=0:0.05:f_range;
fval = f(f_eval);

% Generate sample of path
tic
Xt=(C/b*mR)^(1/a)*GZT(G,Z,Gam,a,t_eval);
toc

%% Compute periodogram and find peaks
Xt_iter=Xt;
Z_hat=[];

MPP=1e-6;
window_name='rect';

tic
r0=[1,1];
for k=1:N
    % zero-padding
%     if length(Xt)>2^13
%         [p,w]=periodogram(Xt_iter);
%     else
%         [p,w]=periodogram(Xt_iter,[],2^13);
%     end
    [p,w]=periodogram(Xt_iter,[],2^15);
    
    [peaks,loc,width]=findpeaks(p,w/delta,'MinPeakProminence',MPP,'SortStr','descend');
    findpeaks(p,w/delta,'MinPeakProminence',MPP)
    
    if isempty(loc) 
        break 
    end
    
    Z_hat(k)=loc(1);
    fun=@(r0,t) r0(1)*cos(t*Z_hat(k))+r0(2)*sin(t*Z_hat(k));
    r=lsqcurvefit(fun,r0,t_eval,Xt_iter',[ ],[ ],optimset('Display','off'));
    Xt_iter=Xt_iter-fun(r,t_eval)';
end
toc

%% save and compute KDE in R
a_str = '150';
n_str = '_0';
% save(['f' num2str(i) 'a' a_str n_str])

% end  % end of loop over all examples 

%% Histogram
nbins=15;
Thist=f_range;
edges=linspace(0,Thist,nbins+1);
figure
histogram(Z_hat,edges,'Normalization','probability')
% histogram([loc;loc2;loc3],nbins,'Normalization','probability')
hold on
histogram(abs(Z),edges,'Normalization','probability','FaceAlpha', 0.4, 'EdgeAlpha',0.4)
hold off
legend('$\hat{Z}$','$\vert Z\vert$','Interpreter','Latex')
title(['Hist. $f_' num2str(i) '$ ($\alpha=' num2str(a) '$, $T= \ ' num2str(T) ...
    '$ , $n= \ ' num2str(length(Xt)') ...
    '$, $N= \ ' num2str(N) '$, $MPP=' num2str(MPP) '$)'], 'interpreter', 'latex','fontsize',10)
