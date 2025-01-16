%-------------------------------------------------------------------------------
%  [function]: simulate states with bilinear interpolation
%-------------------------------------------------------------------------------
function [kt,zt]= simulate_bilinear( kg,zg,kpf, k0,et,phi_z, T0,T1 )

    % set capital policy 2D-interpolant
    [kk,zz]=    ndgrid(kg,zg);
    Kpf=        griddedInterpolant(kk,zz,kpf);      % linear interpolant
    
    % preallocate
    zt= zeros(T1+1,1);
    kt= zeros(T1+1,1);
    % simulate state-paths
    kt(1)= k0;
    for t= 1:T1                                     % loop over time
        % future productivity
        zt(t+1)=    phi_z*zt(t) + et(t);
        % future capital: interpolate policy
        kt(t+1)=    Kpf( kt(t), zt(t) );
    end
    
    % burn initial draws
    zt= zt(T0+1:T1);
    kt= kt(T0+1:T1);
end
