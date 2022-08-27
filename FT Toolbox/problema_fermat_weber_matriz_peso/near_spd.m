function Z=near_spd(X)
    %Clacula una matriz positivia definida cercana a X
    %Sintáxis: Z=near_spd(X)
    Y=0.5*(X+X');
    [Q,D] = eig(Y);
    [n,~]=size(D);
    for i=1:n
        D(i,i)=max([D(i,i) 0]);
    end
    Z=Q*D*Q'+0.1*eye(n);

end
