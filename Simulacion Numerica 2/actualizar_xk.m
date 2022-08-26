function xkm1=actualizar_xk(xk,B,w,s)
  M=entradas_iguales(xk,B);
  [m,p]=size(B);
  xkm1=zeros(m,1);

  if isempty(M)
    vec1=[];
    vec2=1:m;
  else
    vec1=(M(:,1))';
    vec2=setdiff(1:m,vec1);
  end


  for r=1:length(vec1)
    j=M(r,1); i=M(r,2);
    xkm1(j)=B(j,i);
  end

  for j=vec2
    num=0; den=0;
    for i=1:p
      U=w(i)*abs(xk(j)-B(j,i))^(s-2)/norm(xk-B(:,i),s)^(s-1);
      num=num+U*B(j,i);
      den=den+U;
    end
    xkm1(j)=num/den;
  end

end

