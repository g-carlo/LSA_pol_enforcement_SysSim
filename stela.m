function [ x, t, S] = stela(A, b, mu, Max_N)

if nargin<4,
    Max_N=20000;
end;

AA = A'*A;
M = size(A,2);
x = zeros(M,1);
c=0;
precision = 1e-6;

% Set the fast soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% gamma opereator
gamma_op = @(x) min(max(x,0),1);



for t=1:Max_N
    
    x_old=x;
 
    r = diag(AA).*x - AA*x + A'*b;
    Bxt = fast_sthresh(r+c*x,mu)./(diag(AA)+c);
 
    g_p = (A*(Bxt-x))'*(A*(Bxt-x));
    gamma = gamma_op(- ( (A*x-b)'*A*(Bxt-x) + mu*(norm(Bxt,1)-norm(x,1)) )/g_p);

    x = x + gamma*(Bxt-x);

    S = sort(find(x));
    
    err = norm(x_old-x);
    if err< precision
        break;
    end
    
end


end