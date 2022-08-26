function y=grad(B,w,x)
  % Gradiente de la función objetivo en el problema de Fermat-Weber
  % evaluado en x
  % Sintaxis: y=grad(B,w,x)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x = vector de tamaño m
  % Salida:   y = vector de tamaño m, que representa el gradiente evaluado
  %               en x
  % Referencia: - Amir Beck and Shoham Sabach. Weiszfeld’s method: Old and new
  %               results. Journal of Optimization Theory and Applications,
  %               164(1):1–40, 2015.

    m=size(B,2);
    y=0;
    for i=1:m
        y=y+(w(i)/norm(x-B(:,i)))*(x-B(:,i));
    end 
end