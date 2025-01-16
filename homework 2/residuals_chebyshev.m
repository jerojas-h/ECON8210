%-------------------------------------------------------------------------------
%  [function]  Residual function: Chebyshev collocation
%-------------------------------------------------------------------------------
function res= residuals_chebyshev( theta,p, kg, zg,nz,Pi, ...
                                   beta,sigma,eta,alpha,delta )
    % unpack coefficients
    theta=  reshape(theta,p,nz);

    % consumption {t}
    Chat=   chebyshev_poly(kg,p,p)*theta;          % policy at (k,z) nodes
    % get capital policy function 
    [~,Kp]= policy_functions( kg,zg, Chat, sigma,eta,alpha,delta );

    % compute residuals
    res= NaN(p,nz);    
    for iz= 1:nz                                    % loop over z
        % current variables: over (k;z)
        c=      Chat(:,iz);                         % c{t}
        kp=     Kp(:,iz);                           % k{t+1}
        
        % future variables: over (k'(k;z),z')
        % consumption {t+1}
        Cp=     chebyshev_poly( kp, p,p)*theta;
        % return on capital {t+1}
        [~,~,Rp]= policy_functions( kp,zg, Cp, sigma,eta,alpha,delta );

        % EE residuals
        res(:,iz)=  1  -  beta* c.^(sigma) .*(Rp .*Cp.^(-sigma)) *Pi(iz,:)';
    end

    % vectorize for solver
    res= res(:);

end
