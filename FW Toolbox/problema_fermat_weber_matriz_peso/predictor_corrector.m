function [x,k,error]=predictor_corrector(A_array,b_array,x0,y_array,tol)
    % Método Predictor-Corrector para el problema generalizado de Fermat-Weber
    % Sintaxis: predictor_corrector(A_array,b_array,x0,y_array,tol)
    % Entrada:  A_array = [A1^T A2^T ... Ap^T] es una matriz de tamaño n x pm
    %           b_array = [b_1 b_2 ... b_p] es una matriz de tamaño m x p
    %                x0 = vector inicial de tamaño n
    %           y_array = [y_1 y_2 ... y_p] es una matriz de tamaño m x p
    %           tol = constante positiva 
    % Salida:         x = [x_1 x_2 ... x_p] aproximación a la solución del  
    %                     problema generalizado de Fermat Weber
    %             error = error asociado al problema
    %                 k = iteraciones
    % Referencias: - Andersen, K. D., Christiansen, E., Conn, A. R., & Overton, M. L. 
    %                An efficient primal-dual interior-point method for minimizing 
    %                a sum of Euclidean norms. SIAM Journal on Scientific Computing, 
    %                22(1), 243-262. (2000)

    % Paso 1
    x=x0;
    mu=0.1;
    error=tol+1;
    k=0;
    
    while and(error>tol,k<=10000)
        % Paso 2
        k=k+1;
        [m,~,p]=size(A_array);
        z_array=zeros(m,p);
        w_mu=zeros(p,1);
        E_mu=[];
        F_mu=[];
        for i=1:p
            z_array(:,i)=b_array(:,i)-A_array(:,:,i)*x;
            w_mu(i)=sqrt((norm(z_array(:,i)))^2+mu^2);
            E_mu=blkdiag(E_mu,w_mu(i)*eye(m));
            F_mu=blkdiag(F_mu,eye(m)-(1/w_mu(i))*y_array(:,i)*(z_array(:,i))');
        end
        E_mu_inv=pinv(E_mu);
        % Paso 3
        H_mu=0.5*E_mu_inv*(F_mu+F_mu)';
        A=array2block(A_array);
        if rcond(H_mu)<eps
            H_mu=near_spd(H_mu);
        end
        % Paso 4
        z=array2block(z_array);
        y=array2block(y_array);
        delta_x=pinv(A*H_mu*A')*A*E_mu_inv*z;
        delta_z=-A'*delta_x;
        rc=z-E_mu*y;
        delta_y=E_mu_inv*(F_mu*delta_z+rc);        
        delta_y_array=block2array(delta_y,m,p);
        delta_z_array=block2array(delta_z,m,p);
        % Paso 5
        fun_aux=@(beta)fun_beta(z_array,delta_z_array,beta,mu);
        beta = fminbnd(fun_aux,0,1);
        gamma=gamma_max(y_array,delta_y_array);
        x_t=x+beta*delta_x;
        y_t=gamma*(y+delta_y);
        %Si se usa la fórmula alternativa del error, se debe calcular lo siguiente: z_t=z+beta*delta_z;
        % Paso 6
        aux1=block2array(y+delta_y,m,p);
        aux2=block2array(z+delta_z,m,p);
        mu_t=(gap(aux1,aux2))^3/(p*(gap(y_array,z_array))^2);
        w_mu_t=zeros(p,1);
        for i=1:p
            w_mu_t(i)=sqrt((norm(z_array(:,i)))^2+mu_t^2);
        end

        % Paso 7
        rcc=vect_rcc(y_array,z_array,delta_y_array,delta_z_array,w_mu,w_mu_t);
        delta_x=pinv(A*H_mu*A')*A*(E_mu_inv*rcc+y);
        delta_z=-A'*delta_x;
        delta_y=E_mu_inv*(F_mu*delta_z+rcc);

        % Paso 8
        E_mu_t=[];
        for i=1:p
            E_mu_t=blkdiag(E_mu_t,w_mu_t(i)*eye(m));
        end
        vecAux=pinv(E_mu_t)*z;
        if (vecAux)'*delta_z<0
            fun_aux=@(beta)fun_beta(z_array,delta_z_array,beta,mu_t);
            beta = fminbnd(fun_aux,0,1);
            gamma=gamma_max(y_array,delta_y_array);
            x_t=x+beta*delta_x;
            y_t=gamma*(y+delta_y);
            %Si se usa la fórmula alternativa del error, se debe calcular lo siguiente: z_t=z+beta*delta_z;
        end

        % Paso 9
        mu=mu_t;
        xn=x_t;
        y=y_t;
        
        % Paso 10
        y_array=block2array(y,m,p);        
        %Si se usa la fórmula alternativa del error, se debe calcular lo siguiente: 
        %             z=z_t;
        %             z_array=block2array(z,m,p);
        %             error=abs(gap(y_array,z_array));
        error=norm(x-xn);
        x=xn;
    end
end
