function H=hessiano(B,w,x)
  % Hessiano de la función objetivo en el problema de Fermat-Weber
  % evaluado en x
  % Sintaxis: y=hessiano(B,w,x)
  % Entrada:  B = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
  %           w = [w_1 w_2 ... w_p] es un vector de tamaño p
  %           x = vector de tamaño m
  % Salida:   y = matriz de tamaño m x m, que representa el Hessiano evaluado
  %               en x
  % Referencia: - Simone Gorner and Christian Kanzow. On Newton’s method for the
  %               Fermat-Weber location problem. Journal of Optimization Theory
  %               and Applications, 170(1):107–118, 2016.
  
    [d,m]=size(B);
    I=eye(d);
    H=zeros(d);
    for i=1:m
        aux=w(i)/norm(x-B(:,i))^3;
        A=norm(x-B(:,i))^2*I;
        C=(x-B(:,i))*(x-B(:,i))';
        H=H+aux*(A-C);
    end
end