%-------------------------------------------------------------------------------
%  [function]:  guess Chebyshev coefficients
%-------------------------------------------------------------------------------
function theta= guess_chebyshev( kg,p, zg,nz, beta,alpha, opt )

    % guess initial coefficients
    cguess= (1-alpha*beta)*kg.^alpha *exp(zg);
    theta= zeros(p,nz);
    for iz= 1:nz
        % fit policy function guess
        theta(:,iz)= fsolve( @(theta) res_cguess( theta,cguess(:,iz) ), ...
                                                  zeros(p,1), opt );
    end
    % vectorize
    theta= theta(:);

    % [function]  interpolation residuals
    function [res]= res_cguess(theta,cpf)
        % compute residuals
        chat=   chebyshev_poly(kg,p,p)*theta;
        res=    cpf - chat;
    end
end
