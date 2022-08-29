function fun1=fun_beta(z_array,delta_z_array,beta,mu)
    %Calcula la funci�n beta del m�todo predictor-corrector
    %Sintaxis: fun1=fun_beta(z_array,delta_z_array,beta,mu)
    
    [~,p]=size(z_array);
    fun1=0;
    for i=1:p
        fun1=fun1+sqrt(norm(z_array(:,i)+beta*delta_z_array(:,i))^2+mu^2);
    end
end