% Código del Experimento Numérico 3

% Artículo: Aspectos Computacionales del Problema de
%           Fermat-Weber y sus Algoritmos de Solucion

% Autor:  Juan Pablo Soto Quirós
%         Escuela de Matemática
%         Instituto Tecnológico de Costa Rica
%         jusoto@tec.ac.cr

clc; clear;

sigma=2;

%Valores iniciales
%Matrices A1, A2 y A3
A_array(:,:,1)=eye(2); A_array(:,:,2)=sigma*eye(2); A_array(:,:,3)=eye(2);
%Vectores b1, b2 y b3
b_array(:,1)=[-1 0]'; b_array(:,2)=[0 sigma]'; b_array(:,3)=[1 0]';
%Parametros iniciales restantes generados aleatoriamente
y_array=ones(2,3); x=ones(2,1); tol=1e-12;


% Metodo Predictor-Corrector
[xk,k,error_k]=predictor_corrector(A_array,b_array,x,y_array,tol)

%Error Absoluto
error_abs=norm(xk-[0 1]')





