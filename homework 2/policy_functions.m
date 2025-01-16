%-------------------------------------------------------------------------------
%  [function]:  evaluate policy functions with c(k,z)
%-------------------------------------------------------------------------------
function [L,Kp,R]= policy_functions( kg,zg, Chat,  sigma,eta,alpha,delta )

    % labor {t}
    L=      ( (1-alpha)*kg.^alpha *exp(zg) ) .*Chat.^(-sigma);
    L=      L.^( 1/(eta+alpha) );
    % capital{t+1}
    Kp=     kg.^alpha .*L.^(1-alpha) .*exp(zg)  +  (1-delta)*kg  -  Chat;
    % return on capital {t}
    R=      1 + alpha*(kg./L).^(alpha-1) .*exp(zg) - delta;

end
