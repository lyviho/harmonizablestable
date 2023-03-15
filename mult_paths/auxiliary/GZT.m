function Gt=GZT(G,Z,TArr,a,t_eval)
    subsize=1e4;
    subs=max(ceil((length(t_eval)-1)/subsize),1);
    Gt=[];
    if subs>1
       for i=1:subs
         tmp_eval=t_eval((i-1)*subsize+1:i*subsize);
         if i==subs
             tmp_eval=t_eval((i-1)*subsize+1:end);
         end
         [tm,Zm]=ndgrid(tmp_eval,Z);
         Gt=[Gt; (repmat(G(1,:),length(tmp_eval),1).*cos(tm.*Zm)+repmat(G(2,:),length(tmp_eval),1).*sin(tm.*Zm))*(TArr.^(-1/a))'];  
       end
    else
        [tm,Zm]=ndgrid(t_eval,Z);
        Gt=(repmat(G(1,:),length(t_eval),1).*cos(tm.*Zm)+repmat(G(2,:),length(t_eval),1).*sin(tm.*Zm))*(TArr.^(-1/a))';
    end
    
    % Grid for vector evaluation
%     [tm,Zm]=ndgrid(t_eval,Z);
%     Gt=(repmat(G(1,:),length(t_eval),1).*cos(tm.*Zm)+repmat(G(2,:),length(t_eval),1).*sin(tm.*Zm))*(TArr.^(-1/a))';

end