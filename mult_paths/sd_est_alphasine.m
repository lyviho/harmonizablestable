clear 
close all 
clc

%% Examples
% Example 1: N(0,1) pdf for f
stddev=1;
f1=@(x) normpdf(x,0,stddev);
% mR=integral(f1,-inf,inf); 
mR=1;

% Example 2
f2=@(x) 2/sqrt(pi)*x.^2.*exp(-x.^2);
% mR=integral(f,-inf,inf); 
mR=1;
g=@(x) normpdf(x,0,1);
grnd=@() normrnd(0,1);
c_accrej=4*sqrt(2)/exp(1);

%% Series representation
% Set parameters
a=1.5;      % index of stability
b=2.^(a/2).*gamma(1+a./2);
C=(1-a)./(gamma(2-a).*cos(pi*a./2)).*(a~=1)+2/pi*(a==1);

K=1e4;      % number of summands in series representation of X
L=100;      % number of independent paths
n=101;      % number of sample points t_1,...,t_n
T=10;       % sample range [0,R]
t_eval=linspace(0,T,n);

%% Choose example 1 or 2
% ------ % 
% Set ex=1 or ex=2 to choose one of the two examples above. 
% Example 2 takes slightly longer to compute since we use acceptance-
% rejection to generate i.i.d. random variables with density f2
ex=1;       
% ------ %

paths=zeros(L,n);
rng('default')  % set seed for reproducability
for i=1:L
    G=normrnd(0,1,2,K);
    Gam=cumsum(exprnd(1,1,K)); 
    if ex==1
        Z=normrnd(0,stddev,1,K);
    elseif ex==2
        Z=accrejrnd(f2,g,grnd,c_accrej,1,K);
    end
    paths(i,:)=(C/b*mR)^(1/a)*GZT(G,Z,Gam,a,t_eval);
end

figure()
plot(t_eval,paths')
xlabel('$t$','interpreter','latex')
ylabel('$X(t)$','interpreter','latex')
title(['Example ' num2str(ex) ' Paths: $\alpha=$' num2str(a) ', $L=$' num2str(L) ...
        ', $T=$' num2str(T) ', $n=$' num2str(n) ', $K=$' num2str(K)],'interpreter','latex')

%% Koutrouvelis regression type parameter estimation, L indep. paths
% [alpha,beta,sigma,mu]=stablereg(...)
E=zeros(4,n);   % Parameter estimate matrix

% Estimate parameters of X_0
X0=paths(:,1);
[a0,b0,s0,m0] = stablereg(X0);

ca=@(a) integral(@(x) abs(cos(x)).^a,0,2*pi,'ArrayValued',true)/(2*pi);

mR_hat=s0^a0/ca(a0);

% Estimate parameters of X_t-X_0
Xt_X0=paths-repmat(X0,1,n);

for i=2:n
    [E(1,i),E(2,i),E(3,i),E(4,i)] = stablereg(Xt_X0(:,i));
end
at=E(1,:);
st=E(3,:);
 
%% Estimate of codifference function
codiff=2*s0^a0-st.^at;
codiff(1)=codiff(1)+1;    % correcting subtraction of 0^0=1.

%% alpha-sine transform 
if ex==1
    f=f1;
elseif ex==2
    f=f2;
end

Tf=alphasine(f,a);
t2_eval=t_eval/2;
Tf_est=(2*s0^a0-codiff)./(2.^(at+1).*ca(at));

figure()
plot(t2_eval,Tf_est,'Linewidth',2)
hold on
plot(t2_eval,Tf(t2_eval),'--','Linewidth',2)
hold off
title(['$\alpha=$' num2str(a) ', $L=$' num2str(L) ...
    ', $T=$' num2str(T) ', $n=$' num2str(n) ', K=' num2str(K)],'interpreter','latex','fontsize',14)
xlabel('$y$','interpreter','latex','fontsize',14)
ylabel('$Tf(y)$','interpreter','latex','fontsize',14)
legend('$\widehat{Tf}$','$Tf$','interpreter','latex','fontsize',14)

% savefig('ex1Tf.fig')



%% Inversion
[f_est,Ff_est]=aSineInv(Tf_est,T,a0,s0);

figure()
% subplot(1,3,3)
plot(t_eval,f_est(t_eval),'Linewidth',2)
hold on
plot(t_eval,f(t_eval),'--','Linewidth',2)
hold off
% title(['$\alpha=$' num2str(a(j)) ', $N=$' num2str(N) ...
%     ', $R=$' num2str(R) ', $k=$' num2str(K)],'interpreter','latex','fontsize',14)
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$f(x)$','interpreter','latex','fontsize',14)
legend('$\widehat{f}$','$f$','interpreter','latex','fontsize',14)

sgtitle(['$\alpha=$' num2str(a) ', $L=$' num2str(L) ...
        ', $T=$' num2str(T) ', $n=$' num2str(n) ', K=' num2str(K)],'interpreter','latex')


%% Smoothing: Choose mollifier
% Mollifier 1
% e_0=@(x) (1-abs(x)).*(-1<x & x<1);
% e_gamma=@(x,gamma) 1/gamma.*e_0(x./gamma);
psi=@(y) 2*(1-cos(y))./y.^2;

% Mollifier 2
% e_0=@(x) exp(-x.^2)/sqrt(2*pi);
% psi=@(y) exp(-y.^2./(4*pi));

% Set smoothing parameter gamma 
gamma=1;

%% Smoothed Inversion
[f_est,Ff_est]=aSineInvSmooth(Tf_est,T,a0,s0,psi,gamma);

figure()
plot(t_eval,f_est(t_eval),'Linewidth',2)
hold on
plot(t_eval,f(t_eval),'--','Linewidth',2)
hold off
% title(['$\alpha=$' num2str(a(j)) ', $N=$' num2str(N) ...
%     ', $R=$' num2str(R) ', $k=$' num2str(K)],'interpreter','latex','fontsize',14)
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$f(x)$','interpreter','latex','fontsize',14)
legend('$\widehat{f}_{\psi,\gamma}$','$f$','interpreter','latex','fontsize',14)

sgtitle(['$\alpha=$' num2str(a) ', $\gamma=$' num2str(gamma) ', $L=$' num2str(L) ...
        ', $T=$' num2str(T) ', $n=$' num2str(n) ', K=' num2str(K)],'interpreter','latex')
