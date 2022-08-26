function y=norm_canberra(p,q)
    y=sum(abs(p-q)./(abs(p)+abs(q)));
end