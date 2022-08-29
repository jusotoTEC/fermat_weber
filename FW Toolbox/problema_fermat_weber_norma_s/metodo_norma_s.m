function [x,er,k]=metodo_norma_s(B,w,x0,s,tol)
  % Método para el Problema de Fermat-Weber con norma s
  % Sintaxis: [x,error,k]=metodo_norma_s(B,w,x0,s,tol)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x0 = vector inicial de tamaño m
  %           s  = valor real en el intervalo [1,2]
  %           tol = constante positiva 
  % Salida:      x = [x_1 x_2 ... x_p] aproximación a la solución del problema 
  %                  de Fermat Weber
  %          error = error asociado al problema, donde ||x^(k+1)-x^(k)||<tol
  %              k = iteraciones
  % Referencias: - Jack Brimberg and Robert F Love. Global convergence
  %                of a generalized iterative procedure for the minimum
  %                location problem with lp distances. Operations Research,
  %                41(6):1153–1163, 1993

  if or(s>2,s<1)
      error('El valor del parametro s debe estar en el intervalo [1,2]')
  end
  x=x0;
  er=tol+1;
  k=0;
  while and(er>tol,k<=10000)
    x_n=actualizar_xk(x,B,w,s);
    er=norm(x_n-x,s);
    x=x_n;
    k=k+1;
  end

end


