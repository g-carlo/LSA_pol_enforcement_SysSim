function  [x, it, S,g_vec]  = stela3 ( mu, A, y, x_in, eps,gamma, Max_N)
% STELA ALGORITHM for LASSO (PESAVENTO)

% G is the input matrix approximation
% b is the output vector

if nargin<7,
    Max_N=20000;
end;

% Initialization
dim = size(A,2);
x = zeros(dim,1);
PRECISION =1e-6; % old value is 1e-12

AA=A'*A;
d = diag(AA);
c = 0; % usually 0.1

% Set the fast soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% gamma opereator
gamma_op = @(x) min(max(x,0),1);

% % Lambda cost test 
% Lambda_fun = @(x) norm(A*x-y);


W = @ (x) 1./(abs(x)+ eps).^gamma;
g_vec = zeros(1,Max_N);
for t =1:Max_N
    
    x_old=x;
    
    % Main procedure
    r = d.*x - AA*x + A'*y;  % trans
                                                               
    fr = fast_sthresh(r + c*x,mu*W(x_in));
    Bxt = fr./(d + c); % trans                                                                                                                                          
    dif = Bxt -  x;    
                                                                             
    par= dif'*AA*dif; % trans
    gf = -( (x'*A'-y')*A*dif + mu*(norm(W(x_in).*Bxt,1)-norm(W(x_in).*x,1)))/par;
    gamma  = gamma_op(gf); % trans
    g_vec(t) = gamma;
    x = x + gamma*dif;
    
    %cost = Lambda_fun(x);
    
    error = norm(x_old-x);
    it = t;
    
    S = find(x);
        % check if stopping criterion is satisfied or not
    if error < PRECISION 
        break; 
    end;
    
end


end
