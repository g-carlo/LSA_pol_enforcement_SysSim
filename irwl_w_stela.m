function [ x_hat, i, S2, g_vec,error ] = irwl_w_stela(  mu, A, y, x_stela,  eps, gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

PRECISION =1e-6;
Max_N2 = 20;

for i=1:Max_N2
   
    x_old = x_stela;
 
    [x_hat, ~, S2,g_vec] = stela3(mu, A, y, x_stela, eps, gamma);
    
    error = norm(x_old - x_hat);
    x_stela = x_hat;
    if error<PRECISION
        return;
    end
    
end


end

