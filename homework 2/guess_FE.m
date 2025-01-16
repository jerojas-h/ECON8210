%-------------------------------------------------------------------------------
%  [function] FE coefficients guess
%-------------------------------------------------------------------------------
function [theta]= guess_FE( kg,nk, zg,nz,Pi, beta,sigma,eta,alpha,delta, ...
                            do_refined_guess,opt )

    % initial guess: deterministic model without labor
    c0=     (1-alpha*beta)*kg.^alpha *exp(zg);
    theta=  c0(:);

    % refine guess: collocation with B1-splines
    if (do_refined_guess==1)
        th= fsolve( @(theta) err_B1_collocation( theta, kg,nk, zg,nz,Pi), theta,opt );
        theta= reshape(th,nk,nz);
    end
    
    
    % [function]  collocation residuals
    function res= err_B1_collocation( theta, kg,nk, zg,nz,Pi )
        % unpack coefficients
        theta=  reshape(theta,nk,nz);
    
        % consumption {t}
        Chat=   theta;                                  % over (k,z)
        % get capital policy function
        [~,Kp]= policy_functions( kg,zg, Chat, sigma,eta,alpha,delta );

        %  compute residuals
        res= ones(nk,nz);
        % loop over productivity
        for iz= 1:nz
            % current variables: over (k;z)
            c=      Chat(:,iz);                     % c{t}
            kp=     Kp(:,iz);                       % k{t+1}
            
            % future variables: over (k'(k;z),z')
            % consumption {t+1}
            Cp=     interp1( kg,Chat, kp, 'linear','extrap');
            % return on capital {t+1}
            [~,~,Rp]= policy_functions( kp,zg, Cp, sigma,eta,alpha,delta );

            % EE residuals
            res(:,iz)=  1  -  beta* c.^(sigma) .*(Rp .*Cp.^(-sigma)) *Pi(iz,:)';
        end
        % vectorize residuals
        res= res(:);
    end
end
