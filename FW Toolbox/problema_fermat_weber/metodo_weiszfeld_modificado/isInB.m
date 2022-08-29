function [cond,k]=isInB(B,x)
  % Función booleana que indica si el vector x es una columna o no de la
  % matrix B
  % Sintaxis: [cond,k]=isInB(B,x)
  % Entrada:   B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %            x =  es un vector de tamaño m
  % Salida: cond = Condición booleana. cond=True si x es una columna de B.
  %                En caso contrario, entonces cond=False
  %            k = entero no negativo. Si cond=True, entonces j es el
  %                primer número de columna que sea igual a x. Si cond=False,
  %                entonces j=0
  % Referencia: - Amir Beck and Shoham Sabach. Weiszfeld’s method: Old and new
  %               results. Journal of Optimization Theory and Applications,
  %               164(1):1–40, 2015.
  
  p=size(B,2);
  A=x-B;
  errors=zeros(p,1);
  for j=1:p
    errors(j)=norm(A(:,j));
  end
  [min_error,k]=min(errors);
  if min_error<2*eps
    cond=true;
  else
    cond=false;
    k=0;
  end
end

