%-------------------------------------------------------------------------------
% [Function]:  Marginal propensities to consume
%-------------------------------------------------------------------------------
%   Purpose:  computes the model's MPCs
% 
%   [notes]:
%     - the MPC is interpolated for states at which the derivative switches
%     from forward to backward differences
%-------------------------------------------------------------------------------
function [dcpf,cMPC]= MPC( ag,cpf, da,na,ny, t1,nt,A_ay, do_correction,iiF,ii0 )
    
    %---------------------------------------
    %  Instantaneous MPC
    %---------------------------------------
    % compute forward difference
    dcpf=   ( cpf(2:na,:) - cpf(1:na-1,:) )/da;

    %---------------------------------------
    %  MPC over time interval
    %---------------------------------------
    % compute step size
    dt=     t1/(nt-1);

    % initialize
    Ccheck= zeros(na*ny,nt);
    B=      speye(na*ny)/dt - A_ay;
    c=      cpf(:);                                     % vectorize consumption

    % time recursion on Feynman-Kac formula
    for t= nt-1:-1:1
        Ccheck(:,t)=    B\( c + Ccheck(:,t+1)/dt );
    end
    % compute derivate at t=0
    C=      reshape( Ccheck(:,1), na,ny );
    cMPC=   ( C(2:end,:) - C(1:end-1,:) )/da;
    
    %---------------------------------------
    %  interpolate MPCs
    %---------------------------------------
    % patch MPC errors due to switching from forward to backward differences
    if (do_correction==1)

        % loop over productivity
        for iy= 1:ny
           if ( sum(iiF(:,iy))>0 && sum(iiF(:,iy))<na ) % identify troublesome states
        
                % get indexes
                i0=     sum(iiF(:,iy));
                i1=     i0 + sum(ii0(:,iy)) + 1;
                % interpolate MPCs
                ii_int= [i0-1 i1+1];
                dcpf(i0:i1,iy)= interp1( ag(ii_int), dcpf(ii_int,iy), ag(i0:i1));
                cMPC(i0:i1,iy)= interp1( ag(ii_int), cMPC(ii_int,iy), ag(i0:i1));
           end
        end

    end
end