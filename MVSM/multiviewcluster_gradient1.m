function [f,U,iter,F,R,R_0,XXXX] = multiviewcluster_gradient1(B, idx, cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%     B: n-by-n PSD matrix; if it is indefinite, first B=B+aI to make it PSD 
%     idx: an integer vector describing the partition
%     cluster: # of the columns of U
%     U: initial guess;
% Output:
%     f: objective value
%     U: the computed solution
%     r_err_f: relative error of the objective value
%     r_err_g: relative error of the gradient
%     iter: # of iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Incomplete data input!');
elseif nargin < 3
    k = length(idx);
    for i = 1:k
    [Y,~]= qr(randn(idx(i),cluster),0);
    U(sum(idx(1:i-1))+1:sum(idx(1:i)),:)= Y;
    end
    display('the default (5) # of columns is set.');
elseif nargin < 4
    k = length(idx);
    for i = 1:k
    [Y,~]= qr(randn(idx(i),cluster),0);
    U(sum(idx(1:i-1))+1:sum(idx(1:i)),:)= Y;
    end
end

r_err_f = 1;  % (f_k+1 -f_k)/abs(f_k)
r_err_g = 1;  % norm(gradient)
tol_f = 1e-8; tol_g = 1e-4;
maxiter = 1000;
iter = 1;
F = [0];  %f_value
R = [];   
R_0 = []; 
ALPHA = [];  
G = zeros(size(U));
XXXX = [];   
BU = B*U;
value = trace(U'*BU);
F = [F,value];

% -gradient
for i = 1 : k
        Temp = U(sum(idx(1:i-1))+1:sum(idx(1:i)),:)'*BU(sum(idx(1:i-1))+1:sum(idx(1:i)),:);
        Temp = (Temp+Temp')/2;
        G(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = -BU(sum(idx(1:i-1))+1:sum(idx(1:i)),:)...
            +U(sum(idx(1:i-1))+1:sum(idx(1:i)),:)*Temp;
end
lambda = 0.01;
T = U+lambda*G;
R = [R,norm(G)];
% retraction
for i = 1 : k        
        [W,~,V] = svd(T(sum(idx(1:i-1))+1:sum(idx(1:i)),:),0);
        U1(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = W*V';       
end
% select step size of first step
xn = 0;
BU1 = B*U1;
while  value-trace(U1'*BU1) < 0.1*lambda*trace(G'*G) && xn < 100
        xn = xn+1;
        lambda = 0.5*lambda;
        T = U+lambda*G;
        for i = 1 : k
        [W,~,V] = svd(T(sum(idx(1:i-1))+1:sum(idx(1:i)),:),0);
        U1(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = W*V';
        end
end   
BU1 = B*U1;
F = [F,trace(U1'*BU1)];
%% iterating U
while ((r_err_f > tol_f) || (r_err_g > tol_g)) && (iter < maxiter) 
    BU1 = B*U1;
    % the -gradient
    for i = 1 : k
        Temp = U1(sum(idx(1:i-1))+1:sum(idx(1:i)),:)'*BU1(sum(idx(1:i-1))+1:sum(idx(1:i)),:);
        Temp = (Temp+Temp')/2;
        G1(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = -BU1(sum(idx(1:i-1))+1:sum(idx(1:i)),:)...
            +U1(sum(idx(1:i-1))+1:sum(idx(1:i)),:)*Temp;
    end
    r_err_g=norm(G1);
    a = trace((G-G1)'*(U1-U));
    b = trace((G-G1)'*(G-G1));
    alpha = abs(a)/b;
    T = U1+alpha*G1;
    
    % retraction
    for i = 1 : k       
        [W,~,V] = svd(T(sum(idx(1:i-1))+1:sum(idx(1:i)),:),0);
        U2(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = W*V';      
    end
    
    % select step size alpha
    xm = 0;
    while (trace(U1'*(B*U1))-trace(U2'*(B*U2)) < 0.1*alpha*trace(G1'*G1)) && (xm<100)
        xm = xm+1;
        alpha = (1/2)*alpha;
        T = U1+alpha*G1;
        for i = 1 : k
        [W,~,V] = svd(T(sum(idx(1:i-1))+1:sum(idx(1:i)),:),0);
        U2(sum(idx(1:i-1))+1:sum(idx(1:i)),:) = W*V';
        end
    end
    U = U1;
    G = G1;
    U1 = U2;
    ALPHA = [ALPHA,alpha];
    f = trace(U1'*BU1);
    r_err_f=(f-F(end))/abs(f);
    F = [F, f];
    R_0 = [R_0,r_err_f];
    R = [R, r_err_g];
    XXXX=[XXXX,xm];
    iter = iter+1;
end

