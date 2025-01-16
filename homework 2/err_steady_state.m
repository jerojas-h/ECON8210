%-------------------------------------------------------------------------------
%  [function]  steady state errors
%-------------------------------------------------------------------------------
function err= err_steady_state( x, beta,sigma,psi,alpha,delta )

    % unpack variables
    c= x(1);  l= x(2);  k= x(3); 

    % preallocate
    err= ones(3,1);
    % euler equation
    err(1)=     1 - beta*( alpha*(k/l)^(alpha-1) + (1-delta) );
    % intratemporal foc
    err(2)=     l^psi - (1-alpha)*(k/l)^alpha *c^(-sigma);
    % budget constraint
    err(3)=     c + k -  k^alpha*l^(1-alpha) - (1-delta)*k;

end

