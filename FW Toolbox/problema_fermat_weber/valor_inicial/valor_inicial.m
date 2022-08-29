function x0=valor_inicial(B,w)
  % Valor inicial de los métodos iterativos que aproxima la solución del
  % problema de Fermat-Weber
  % Sintaxis: x0=valor_inicial(B,w)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  % Salida:  x0 = vector de tamaño m
  % Referencia: - Amir Beck and Shoham Sabach. Weiszfeld’s method: Old and new
  %               results. Journal of Optimization Theory and Applications,
  %               164(1):1–40, 2015.
    p=size(B,2);
    f_b=zeros(p,1);
    for i=1:p
        f_b(i)=funcion_objetivo(B,w,B(:,i));
    end
    [~,j]=min(f_b);
    Rj=R_j(B,w,j);
    dj=(-1/norm(Rj))*Rj;
    Lj=L_j(B,w,j);
    alpha_j=(norm(Rj)-w(j))/Lj;
    x0=B(:,j)+alpha_j*dj;
end
