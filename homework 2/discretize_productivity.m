%-------------------------------------------------------------------------------
%  [function]:  AR(1) discretization
%-------------------------------------------------------------------------------
%   Purpose:    discretizes a de-meaned AR(1) process using Tauchen's
%               method
%-------------------------------------------------------------------------------
function [zg,Pi]= discretize_productivity( nz, phi,sde, m )

    % build productivity grid
    sdz=    sde/sqrt(1-phi^2);                      % productivity std 
    zmin=   -m*sdz;                                 % grid lower bound
    zmax=   m*sdz;                                  % grid upper bound
    dz=     (zmax-zmin)/(nz-1);                     % step size
    zg=     zmin + ((1:nz)-1)*dz;                   % productivity grid
    
    % build transition matrix
    Pi= normcdf( (zg+dz/2 - phi*zg')/sde ) - normcdf( (zg-dz/2 - phi*zg')/sde );
    % adjust tails 
    Pi(:,1)=    normcdf( (zg(1)+dz/2 - phi*zg')/sde );
    Pi(:,nz)=   1 - normcdf( (zg(nz)-dz/2 - phi*zg')/sde );

end
