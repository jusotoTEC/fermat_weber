% C�digo del Experimento Num�rico 2

% Art�culo: Aspectos Computacionales del Problema de
%           Fermat-Weber y sus Algoritmos de Solucion

% Autor:  Juan Pablo Soto Quir�s
%         Escuela de Matem�tica
%         Instituto Tecnol�gico de Costa Rica
%         jusoto@tec.ac.cr

clc; clear;

tam=[10 50 100 200 400 600 800 1000];
vec_s=[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];

time=zeros(6,length(tam));
iter=zeros(6,length(tam));
error_abs=zeros(6,length(tam));

m=20;
tol=eps;


cont=0;
for p=tam
    display(['Ejecutando p = ' num2str(p)])
    cont=cont+1;
    B=randn(m,p);
    w=rand(p,1);

    % Valor inicial
    x0=valor_inicial(B,w);
   
    %M�todo de Norma s
    for r=1:length(vec_s)
        tic; [x,error,k]=metodo_norma_s(B,w,x0,vec_s(r),tol); t=toc;
        time(r,cont)=t; iter(r,cont)=k; error_abs(r,cont)=error;
    end
end

%Gr�ficas

%Dimensi�n vrs Errores

figure 
hold on
for r=1:length(vec_s)
    plot(tam,error_abs(r,:),'DisplayName',['s=', num2str(vec_s(r))])
end
grid on
lgd = legend;
lgd.NumColumns = 2;
title('Dimensi�n vrs Error')
xlabel('Dimension (p)')
ylabel('Error')

%Dimensi�n vrs Iteraciones

figure 
hold on
for r=1:length(vec_s)
    plot(tam,iter(r,:),'DisplayName',['s=', num2str(vec_s(r))])
end
grid on
lgd = legend;
lgd.NumColumns = 2;
title('Dimensi�n vrs Iteraciones')
xlabel('Dimension (p)')
ylabel('Iteraciones')

%Dimensi�n vrs Tiempo

figure 
hold on
for r=1:length(vec_s)
    plot(tam,time(r,:),'DisplayName',['s=', num2str(vec_s(r))])
end
grid on
lgd = legend;
lgd.NumColumns = 2;
title('Dimensi�n vrs Tiempo')
xlabel('Dimension (p)')
ylabel('Tiempo')


