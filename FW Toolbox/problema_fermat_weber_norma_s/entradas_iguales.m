function M=entradas_iguales(x,B)
  A=x-B;
  B=abs(A);
  [j_v,i_v]=find(B<eps);
  if isempty(j_v)
    M=[];
  else
    m=length(j_v);
    M=[j_v(1) i_v(1)];
    for k=2:m
      if ~ismember(j_v(k),M(:,1))
        M=[M; j_v(k) i_v(k)];
      end
    end
  end
end

