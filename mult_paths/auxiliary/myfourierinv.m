function Finvf=myfourierinv(f,lower,upper)
%     kernel=@(x,y) exp(-1i.*x.*y);
    kernel=@(x,y) cos(x.*y);
%     Finvf=@(y) 1/pi*integral(@(x) kernel(x,y).*f(x),0,upper,'AbsTol',1e-3,'RelTol',1e-3,'ArrayValued',true);
%     Finvf=@(y) 1/(2*pi)*integral(@(x) kernel(x,y).*f(x),lower,upper,'AbsTol',1e-3,'RelTol',1e-3,'ArrayValued',true);
    Finvf=@(y) 1/(2*pi)*integral(@(x) kernel(x,y).*f(x),lower,upper,'AbsTol',1e-3,'RelTol',1e-3,'ArrayValued',true);
end