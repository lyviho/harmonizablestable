function c=c_coeff(alpha,j)
%     if sumj==0
%         a=1/2^(alpha/2)*gamma(1/2+alpha/2)/(gamma(1+alpha/4)*gamma(1/2+alpha/4));
%     else
%         a=(-1).^j./2.^(alpha-1).*gamma(1+alpha)./(gamma(1+alpha/2-j).*gamma(1+alpha/2+j));
%     end
% c=(-1).^j./2.^(alpha-1).*gamma(1+alpha)./(gamma(1+alpha/2-j).*gamma(1+alpha/2+j));
c=(-1).^j./2.^(alpha).*gamma(1+alpha)./(gamma(1+alpha/2-j).*gamma(1+alpha/2+j));
end