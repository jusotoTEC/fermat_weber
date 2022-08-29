function [x,error,k]=metodo_weiszfeld(B,w,x0,tol)
  % Método Weiszfeld
  % Sintaxis: [x,error,k]=metodo_weiszfeld(B,w,x0,tol)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x0 = vector inicial de tamaño m
  %           tol = constante positiva 
  % Salida:      x = [x_1 x_2 ... x_p] aproximación a la solución del problema 
  %                  de Fermat Weber
  %          error = error asociado al problema, donde ||x^(k+1)-x^(k)||<tol
  %              k = iteraciones
  % Referencias: - Endre Weiszfeld. Sur le point pour lequel la somme des
  %                distances de n points donnes est minimum. Tohoku Mathematical
  %                Journal, First Series, 43:355–386, 1937
  %              - Simone Gorner and Christian Kanzow. On Newton’s method for
  %                the Fermat-Weber location problem. Journal of Optimization
  %                Theory and Applications, 170(1):107–118, 2016

  p=size(B,2);
  f_b=zeros(p,1);
  for i=1:p
      f_b(i)=funcion_objetivo(B,w,B(:,i));
  end
  [~,j]=min(f_b);
  nRj=norm(R_j(B,w,j));
  if nRj<=w(j)
    x=B(:,j);
    error=0;
  else
    error=tol+1;
    k=0;
    x=x0;
    while and(error>tol,k<=10000)
      alpha=0;
      y=0;
      for i=1:p
          aux=w(i)/norm(x-B(:,i));
          alpha=alpha+aux;
          y=y+aux*B(:,i);
      end
      x_n=(1/alpha)*y;
      error=norm(x_n-x);
      x=x_n;
      k=k+1;
    end
  end
end
