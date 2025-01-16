%-------------------------------------------------------------------------------
%  [function]  evaluate tent basis
%-------------------------------------------------------------------------------
function [psi_hat]= tent_basis( xg,nx, kg,n )

    % preallocate
    psi_hat= zeros(nx,n);

    % compute basis elements over x-grid
    for i= 1:n                                      % loop over elements
        if (i==1)
            % get interval bounds
            k2= kg(i+1);  k1= kg(i);
            % build basis
            ii=             (xg<=k2);
            psi_hat(:,i)=   (k2-xg)/(k2-k1).*ii;
        elseif (i<n)
            % get interval bounds
            k2= kg(i+1);  k1= kg(i);  k0= kg(i-1);
            % build basis
            ii1=            (xg>k0).*(xg<=k1);
            ii2=            (xg>k1).*(xg<=k2);   
            psi_hat(:,i)=   (xg-k0)/(k1-k0).*ii1  +  (k2-xg)/(k2-k1).*ii2;
        else
            % get interval bounds
            k1= kg(i);  k0= kg(i-1);
            % build basis
            ii1=            (xg>k0);
            psi_hat(:,i)=   (xg-k0)/(k1-k0).*ii1;
        end
    end
end
