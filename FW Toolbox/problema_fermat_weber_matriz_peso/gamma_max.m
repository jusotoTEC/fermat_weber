function gamma_min=gamma_max(y_array,delta_y_array)    
    %Calcula un valor m�ximo de gama del m�todo predictor-corrector
    %Sint�xis: gamma_min=gamma_max(y_array,delta_y_array)    
    [~,p]=size(y_array);
    gamma=zeros(p,1);
    for i=1:p
        gamma(i)=1/(norm(y_array(:,i)+delta_y_array(:,i)));
    end        
    gamma_min=min(gamma);
end