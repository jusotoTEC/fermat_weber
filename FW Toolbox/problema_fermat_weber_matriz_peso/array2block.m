function Z=array2block(Z_array)
    % Convierte un arreglo de dimension m x n x p en una matriz de tamaño
    % mp x n, concatenando matriz por matriz
    % Sintáxis: Z=array2block(Z_array)


    tam=size(Z_array);
    if length(tam)==2
        [~,p]=size(Z_array);
        Z=[];
        for i=1:p
          Z=[Z; Z_array(:,i)];
        end
    elseif length(tam)==3
        [~,~,p]=size(Z_array);
        Z=[];
        for i=1:p
          Z=[Z (Z_array(:,:,i))'];
        end
    end
end