% Código del Experimento Numérico 4

% Artículo: Aspectos Computacionales del Problema de
%           Fermat-Weber y sus Algoritmos de Solucion

% Autor:  Juan Pablo Soto Quirós
%         Escuela de Matemática
%         Instituto Tecnológico de Costa Rica
%         jusoto@tec.ac.cr

clc; clear; close all


% Paso 1: Cargar imágenes de la base de datos

numBD=1; % Base de rostros con 12 imágenes
%numBD=2; % Base de rostros con 400 imágenes

a=dir(['base de datos de caras ' num2str(numBD) '\' '/*.jpg']);
numImagenes=size(a,1);


[m,n]=size(im2double(imread(['base de datos de caras ' num2str(numBD) '\cara (1).jpg'])));


B=[];
for i=1:numImagenes
    text=['base de datos de caras ' num2str(numBD) '\cara (' num2str(i) ').jpg'];
    Aux=im2double(imread(text));
    B=[B Aux(:)];    
end


% Paso 2: Calcular imagen promedio

imgPromedio=mean(B,2);

% Paso 3: Calcular las distancias

dist=[];
for k=1:numImagenes
    dist=[dist norm(imgPromedio-B(:,k))];
end

figure
bar(1:numImagenes,dist)
xlabel('Número de Imagen (j)')
ylabel('Distancia')

% Paso 4: Calcular los pesos wi

w=1-0.9*dist/max(dist);

figure
bar(1:numImagenes,w)
xlabel('Número de Peso (j)')
ylabel('Peso \omega_j')


%Paso 5: Calcular imagen prototipo usando el método de Weiszfeld

x0=valor_inicial(B,w);
[imgPrototipo,er1,k1]=metodo_weiszfeld(B,w,x0,10e-10);
er2=norm(grad(B,w,imgPrototipo));


figure
imshow(reshape(imgPrototipo,[m,n]))

er_ssim=ssim(reshape(imgPrototipo,[m,n]),reshape(imgPromedio,[m,n]))
er_mse=norm(imgPrototipo-imgPromedio)
