%-------------------------------------------------------------------------------
%  [function]  Euler equation residuals
%-------------------------------------------------------------------------------
function [ee_res]=  euler_eqn_residuals( cpf,kpf, method, zg,nz, Pi, ...
                                         beta,sigma,eta,alpha,delta )

    % evaluate residuals
    ee_res= NaN(size(cpf));
    for iz= 1:nz                                    % loop over z
        % current variables: over (k;z)
        c=      cpf(:,iz);                          % c{t}
        kp=     kpf(:,iz);                          % k{t+1}
        
        % future variables: over (k'(k;z),z')
        % consumption {t+1}
        if (method.method==1)
            Cp= chebyshev_poly( kp,length(kp), method.p )*method.theta;
        elseif (method.method==2)
            Cp= tent_basis( kp,length(kp), method.kg,method.n )*method.theta;
        elseif (method.method==3)
            Cp= perturbed_pf( kp-method.kss,zg, method.DR,method.idx_c, method.order );
        end
        % return on capital {t+1}
        [~,~,Rp]= policy_functions( kp,zg, Cp, sigma,eta,alpha,delta );

        % EE residuals
        inverse_ee=     (  beta *(Rp .*Cp.^(-sigma)) *Pi(iz,:)'  ).^(-1/sigma);
        ee_res(:,iz)=   inverse_ee./c - 1;
    end
end

