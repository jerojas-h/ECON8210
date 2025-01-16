%-------------------------------------------------------------------------------
%  [function]: evaluate perturbed policy function
%-------------------------------------------------------------------------------
%   [notes]: 
%     - constants are in 'declared order' and higher-order terms are in
%       'decision rule' order
%-------------------------------------------------------------------------------
function pf= perturbed_pf( dK,dZ, DR,i, approx_order )

    % get steady-states
    G0= DR.ys(DR.order_var);

    % get stochastic steady-state correction
    if (approx_order>1)
        GS2= DR.ghs2(DR.order_var);
    else
        GS2= 0*G0;
    end

    % get first-order coefficients
    G1= DR.ghx;
    
    
    %------------------------------------------
    % compute policy functions
    %------------------------------------------
    % constants
    constants=  G0(i) + 1/2*GS2(i);

    % 1st-order terms
    g_k=    G1(i,1)*dK;
    g_z=    G1(i,2)*dZ;
    order1= g_k + g_z;
   
    % 2nd-order terms
    order2= 0;
    if (approx_order>1)
        % get coefficients
        G2=     DR.ghxx;
        % build approximation terms
        g_k2=   G2(i,1) *dK.^2;
        g_kz=   1/2*sum( G2(i,[2 3]) ) *dK *dZ;
        g_z2=   G2(i,4) *dZ.^2;
        % 2nd-order terms
        order2= 1/2*(g_k2 + g_kz + g_z2);
    end

    % 3rd-order terms
    order3= 0;
    if (approx_order>2)
        % get coefficients
        G3=     DR.ghxxx;
        % build approximation terms
        g_k3=   G3(i,1) *dK.^3;
        g_zk2=  1/3* sum( G3(i,[2 3 5]) ) *dK.^2 *dZ;
        g_kz2=  1/3* sum( G3(i,[4 6 7]) ) *dK *dZ.^2;
        g_z3=   G3(i,8) *dZ.^3;
        % 3rd order term
        order3= 1/6*(g_k3 + g_zk2 + g_kz2 + g_z3);
    end    
    
    % policy function
    pf= constants + order1 + order2 + order3;

end
