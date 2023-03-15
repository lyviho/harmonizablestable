%----------------------------------------------------------------%
%-- @usage:         T( f , alpha)                              --%
%-- @return:        function handle                            --%
%-- @param:         f - even function                          --%
%--                 alpha - value of alpha in (0,2]            --%
%-- @description:   Integral transform T                       --%
%----------------------------------------------------------------%
function Tf=alphasine(f,alpha)
%     kernel=@(x,t) abs(sin(t.*x/2)).^alpha;
    kernel=@(x,t) abs(sin(t.*x)).^alpha;
    Tf=@(t) integral(@(x) kernel(x,t).*f(x),0,inf,'ArrayValued',true);
end