function g=gap(y_array,z_array)
    %Calcula la función gap
    %Sintaxis: g=gap(y_array,z_array)
    [~,p]=size(z_array);
    g=0;
    for i=1:p
        g=g+(norm(z_array(:,i))-(y_array(:,i))'*z_array(:,i));
    end
end