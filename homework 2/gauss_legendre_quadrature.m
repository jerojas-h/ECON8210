%-------------------------------------------------------------------------------
%  [function]  Gauss-Legendre quadrature
%-------------------------------------------------------------------------------
function [x,w]= gauss_legendre_quadrature(a,b)
    % [note]: computes Guass-Legendre nodes & weights for n=10 nodes over
    % the interval [a,b]. Roots from numerical recipes in C

    % nodes
    x=  [ 0.1488743389  0.4333953941  0.6794095682  0.8650633666  0.9739065285 ];
    x=  [ -x(5:-1:1) x ]';
    % weights
    w=  [ 0.2955242247  0.2692667193  0.2190863625  0.1494513491  0.0666713443 ];
    w=  [ w(5:-1:1) w ]';

    % transform to [a,b] interval
    x=  (b-a)/2 *x  +  (a+b)/2;
    w=  (b-a)/2 *w;
end

