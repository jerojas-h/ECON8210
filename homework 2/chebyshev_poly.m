%-------------------------------------------------------------------------------
%  [function]:  evaluate Chebyshev polynomials
%-------------------------------------------------------------------------------
function T= chebyshev_poly(xg,nx,m)

    % build polynomial
    T=      ones(nx,m);                             % order 0
    T(:,2)= xg;                                     % order 1
    for p= 3:m
        % recursion for order p>1
        T(:,p)= 2*xg.*T(:,p-1) - T(:,p-2);          
    end    
end
