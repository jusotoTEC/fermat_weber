function Lj=L_j(B,w,j)
  % Operador L en la ecuación (12) en el artículo "Weiszfeld’s method: Old 
  % and new results", cuando x=b_j
  % Sintaxis: Lj=L_j(B,w,j)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p  
  %           j = entero positivo entre 1 y p
  % Salida: Lj = constante positiva
  % Referencia: - Amir Beck and Shoham Sabach. Weiszfeld’s method: Old and new
  %               results. Journal of Optimization Theory and Applications,
  %               164(1):1–40, 2015.
  
  p=size(B,2);
  Lj=0;
  for i=1:p
    if i~=j
      Lj=Lj+(w(i)/norm(B(:,j)-B(:,i)));
    end
  end
end

