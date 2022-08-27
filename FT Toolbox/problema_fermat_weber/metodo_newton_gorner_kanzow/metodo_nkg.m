function [x,error,k]=metodo_nkg(B,w,x0,tol)
  % Método Newton-Gorner-Kanzow
  % Sintaxis: [x,error,k]=metodo_nkg(B,w,tol)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x0 = vector inicial de tamaño m
  %           tol = constante positiva 
  % Salida:     x = [x_1 x_2 ... x_p] aproximación a la solución del problema 
  %               de Fermat Weber
  %         error = error asociado al problema, donde ||x^(k+1)-x^(k)||<tol
  %             k = iteraciones
  % Referencia: - Simone Gorner and Christian Kanzow. On Newton’s method for the
  %               Fermat-Weber location problem. Journal of Optimization Theory
  %               and Applications, 170(1):107–118, 2016.

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
    rho=rand(1); 
    sigma=rand(1)/2;
    error=tol+1;
    k=0;
    x=x0;
    while error>tol
        k=k+1;
        M=hessiano(B,w,x);
        z=-grad(B,w,x);
        d=linsolve(M,z);
        alpha=paso(B,w,x,d,rho,sigma);
        x_n=x+alpha*d;
        error=norm(x_n-x);
        x=x_n;
    end
  end
end
