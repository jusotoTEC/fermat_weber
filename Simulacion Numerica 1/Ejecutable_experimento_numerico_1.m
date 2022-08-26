% Código del Experimento Numérico 1

% Artículo: Aspectos Computacionales del Problema de
%           Fermat-Weber y sus Algoritmos de Solucion

% Autor:  Juan Pablo Soto Quirós
%         Escuela de Matemática
%         Instituto Tecnológico de Costa Rica
%         jusoto@tec.ac.cr

clc; clear;

tam=[10 100 500 750 1000 2000 5000 10000 15000 20000 30000];

time=zeros(3,length(tam));
iter=zeros(3,length(tam));
error_abs=zeros(3,length(tam));
error_grad=zeros(3,length(tam));

m=20; 
tol=eps;

k=0;
for p=tam
    display(['Ejecutando p = ' num2str(p)])
    k=k+1;    
    B=randn(m,p);
    w=rand(p,1);
    
    % Valor inicial
    x0=valor_inicial(B,w);
    
    %Método de Weiszfeld
    tic; [x1,error1,k1]=metodo_weiszfeld(B,w,x0,tol); t1=toc;  e1=norm(grad(B,w,x1));    
    time(1,k)=t1; iter(1,k)=k1; error_abs(1,k)=error1; error_grad(1,k)=e1;
    
    %Método de Weiszfeld Modificado
    tic; [x2,error2,k2]=metodo_weiszfeld_mod(B,w,x0,tol); t2=toc; e2=norm(grad(B,w,x2));
    time(2,k)=t2; iter(2,k)=k2; error_abs(2,k)=error2; error_grad(2,k)=e2;
    
    %Método de Newton-Gorner-Kanzow 
    tic; [x3,error3,k3]=metodo_nkg(B,w,x0,tol); t3=toc; e3=norm(grad(B,w,x3));
    time(3,k)=t3; iter(3,k)=k3; error_abs(3,k)=error3; error_grad(3,k)=e3;
end

% Construcción de Tablas con Información Numérica de los Experimentos

  dim_tam=length(tam);
  dimension=cell(1,dim_tam);
  
  for k=1:dim_tam
      dimension{k}=num2str(tam(k));
  end

%Tabla 1 =  Tiempos e Iteraciones
fprintf('Tabla 1: Tiempos e Iteraciones\n')
table_Results=table((time(1,:))',(time(2,:))',(time(3,:))',(iter(1,:))',(iter(2,:))',(iter(3,:))', 'RowNames', dimension);
table_Results.Properties.VariableNames={'Time_Weiszfeld', 'Time_Weiszfeld_Mod', 'Time_NKG','Iteraciones_Weiszfeld', 'Iteraciones_Weiszfeld_Mod', 'Iteraciones_NKG'};
disp(table_Results)  
  
%Tabla 2 =  Errores Absolutos y Gradiente
fprintf('Tabla 2: Errores Absolutos y Gradiente\n')
table_Results=table((error_abs(1,:))',(error_abs(2,:))',(error_abs(3,:))',(error_grad(1,:))',(error_grad(2,:))',(error_grad(3,:))', 'RowNames', dimension);
table_Results.Properties.VariableNames={'Error_Abs_Weiszfeld', 'Error_Abs_Weiszfeld_Mod', 'Error_Abs_NKG','Error_Gra_Weiszfeld', 'Error_Gra_Weiszfeld_Mod', 'Error_Gra_NKG'};
disp(table_Results)
