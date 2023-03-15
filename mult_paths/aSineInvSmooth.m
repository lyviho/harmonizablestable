function [f_est,Ff_est]=aSineInvSmooth(Tf_sample,R,a,s,psi,gamma)

N=length(Tf_sample);
c=c_coeff(a,1:N);
C=zeros(N,N);
for i=1:N
    tmp=[zeros(1,i-1),reshape([c;zeros(i-1,length(c))],1,[])];
    tmp=tmp(1:N);
    C(i,:)=tmp;
end

ca=integral(@(x) abs(cos(x)).^a,0,2*pi)/(2*pi);
mR=s^a/ca;

eta=Tf_sample-c_coeff(a,0)/2*mR;
xi=C\eta';
fhat_sample=[flip(xi);mR;xi];
% f_est=@(x_eval) bandliminv(fhat_sample,N,R,x_eval);   %R/2?!!!
Ff_est=@(bli_eval) bandlim(fhat_sample,N,R,bli_eval);
f_est=myfourierinv(@(y) Ff_est(y).*psi(gamma*y),-R,R);   

end