function rcc=vect_rcc(y_array,z_array,delta_y_array,delta_z_array,w_mu,w_mu_t)
    %Calcula el vector r_c^c del método predictor-corrector
    %Sintáxis: rcc=vect_rcc(y_array,z_array,delta_y_array,delta_z_array,w_mu,w_mu_t)
    [~,p]=size(y_array);
    rcc=[];
    for i=1:p
        p1=z_array(:,i);
        p2=w_mu_t(i)*y_array(:,i);
        p3=(((z_array(:,i))'*delta_z_array(:,i))/(w_mu_t(i)))*delta_y_array(:,i);
        p4=(norm(delta_z_array(:,i))^2/(2*w_mu_t(i)))*y_array(:,i);
        p5=(((z_array(:,i))'*delta_z_array(:,i))^2/(2*(w_mu_t(i))^3))*y_array(:,i);
        p=p1-p2-p3-p4+p5;        
        
        q1=(w_mu(i)-w_mu_t(i))*delta_y_array(:,i);
        q2=(1/w_mu(i)-1/w_mu_t(i))*((z_array(:,i))'*delta_z_array(:,i))*y_array(:,i);
        q=q1+q2;
        
        t=(1/(2*w_mu(i)))*(((delta_z_array(:,i))'*y_array(:,i))*z_array(:,i)-((delta_z_array(:,i))'*z_array(:,i))*y_array(:,i));
        
        rcc=[rcc; p+q+t];
    end
end