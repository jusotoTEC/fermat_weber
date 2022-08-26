function Z_array=block2array(Z,m,p)
    %Convierte una matriz en un vector columna
    %Sintáxis: Z_array=block2array(Z,m,p)
    
    Z_array=zeros(m,p);
    for i=1:p
        Z_array(:,i)=Z((i-1)*m+1:i*m);
    end
end