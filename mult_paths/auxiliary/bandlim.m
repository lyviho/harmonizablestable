function Ff=bandlim(Ff_sol,N,R,bli_eval)
    [x,n] = ndgrid(bli_eval,-N:N);
    Ff = sinc((N/R*x - n))*Ff_sol;
%     Ff=@(y) sinc((y-y_sample)/(N/L))*Ff_sol;
end