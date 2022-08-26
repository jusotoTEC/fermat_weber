function alpha=paso(B,w,x,d,rho,sigma)
  % Tamaño de paso en el Algoritmo Newton-Gorner-Kanzow
  % Sintaxis: alpha=paso(B,w,x,rho,d,sigma)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x = vector de tamaño m
  %           d = vector de tamaño m
  %         rho = variable aleatoria en el intervalo ]0,1[   
  %       sigma = variable aleatoria en el intervalo ]0,0.5[   
  % Salida:  alpha = constante positiva
  % Referencia: - Simone Gorner and Christian Kanzow. On Newton’s method for
  %               the Fermat-Weber location problem. Journal of Optimization
  %               Theory and Applications, 170(1):107–118, 2016
  
    cond=1;
    k=0;
    while cond==1
        alpha=rho^k;
        l=funcion_objetivo(B,w,x+alpha*d);
        r=funcion_objetivo(B,w,x)+sigma*alpha*(grad(B,w,x))'*d;
        if l<=r
            cond=0;
        else
            k=k+1;            
        end
    end        
end