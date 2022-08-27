function f=funcion_objetivo(B,w,x)
  % Funci�n objetivo del problema Fermat-Weber evaluado en x
  % Sintaxis: f=funcion_objetivo(B,w,x)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tama�o m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tama�o p
  %           x = vector de tama�o m
  % Salida:   f = n�mero no negativo. 
  % Referencia: - Amir Beck and Shoham Sabach. Weiszfeld�s method: Old and new
  %               results. Journal of Optimization Theory and Applications,
  %               164(1):1�40, 2015.
  
  p=size(B,2);
  f=0;
  for i=1:p
    f=f+w(i)*norm(x-B(:,i));
  end
end

