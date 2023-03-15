function f=bandliminv(Ff_sol,N,R,inv_eval)
    [x,n] = ndgrid(inv_eval,-N:N);
%     f=(1/(2*pi*R)*rectangularPulse(x/(2*pi*R)).*exp(-1i*n/R.*x))*Ff_sol;
    f=(1/(2*pi.*N./R).*rectangularPulse(x./(2.*pi.*N./R)).*exp(-1i.*n./(N./R).*x))*Ff_sol;
%     f=(1/(pi*N/R)*rectangularPulse(x/(pi*N/R)).*exp(-1i*n/(N/R).*2.*x))*Ff_sol;
end