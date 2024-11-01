clear variables; close all; clc;

% program options
store_results=  0;
load_results=   1;

% parameters
beta=   0.97;
delta=  0.1;
alpha=  0.33;
phi=    2*0.05;

% grid parameters
nK=     7;
nI=     7;
nk=     51;
ni=     51;
k_fac=  0.3;
i_fac=  0.5;


% value function parameters
tol_vf=             1e-4;
max_iters_vf=       400;
% price schedule parameters
tol_LoM=            1e-2;
max_iters_LoM=      40 + 1;

% specs
do_deterministic=   0;


% tax grid
tg= [0.2 0.25 0.3]';
nt= length(tg);

% transition matrices
Pi_tau= [ 0.9   0.1  0
          0.05  0.9  0.05
          0     0.1  0.9  ];

% productivity grid
zg= exp([-0.0673  -0.0336 0   0.0336  0.0673]');
nz= length(zg);
Pi_z= [ 0.9727  0.0273  0       0       0
        0.0041  0.9806  0.0153  0       0
        0       0.0082  0.9836  0.0082  0
        0       0       0.0153  0.9806  0.0041
        0       0       0       0.0273  0.9727 ];

% combined transition: (z,tau)
Pi= kron(Pi_z,Pi_tau);
ns= nt*nz;



%----------------------------------------------------------------------------
%  steady states
%----------------------------------------------------------------------------
% parameters
tss=    0.25;
zss=    0;

% discount rate
rho=    1/beta-1;

% interest rate

% labor, capital, & investment
lss=    sqrt(  1/( 1 + rho/((rho+delta)*(1-tss)) * alpha/(1-alpha) )  );
kss=    ( alpha/(rho+delta)*exp(zss) )^(1/(1-alpha)) *lss;
iss=    delta*kss;

% interest rate & wage
rss=    rho+delta;
wss=    (1-alpha)*exp(zss)* (kss/lss)^alpha;

% consumption
css=    (1-tss)*wss/lss;
gss=    tss*wss*lss;
yss=    exp(zss)* kss^alpha *lss^(1-alpha);

% grids
kg= linspace( (1-k_fac)*kss, (1+k_fac)*kss, nk)';
ig= linspace( (1-i_fac)*iss, (1+i_fac)*iss, ni)';

% aggregate grids
Kg_agg= linspace( (1-k_fac)*kss, (1+k_fac)*kss, nK)';
Ig_agg= linspace( (1-i_fac)*iss, (1+i_fac)*iss, nI)';



%----------------------------------------------------------------------------
%  solve the recursive competitive equilibrium
%----------------------------------------------------------------------------
% do model under certainty
if (do_deterministic==1)
    tg= tss;  nt= 1;
    zg= 1;    nz= 1;
    Pi= 1;    ns=  1;
end


% investment adjustment-cost over (i-,i)
Phi= (  1 - phi/2*( ig'./ig -1 ).^2  ).*ig';

% laws of motion
Hk0= ones(ns,nK,nI)*kss;
Hi0= ones(ns,nK,nI)*iss;
Hl0= ones(ns,nK,nI)*lss;   


% value function guess
VV0= zeros(ni*ns,nk,nK,nI);
[xx,yy,zz,ww]= ndgrid(ig,kg,Kg_agg,Ig_agg);

% load results
if (load_results==1)
    load('Hk0.mat');
    load('Hi0.mat');
    load('Hl0.mat');
    load('VV0.mat');
end



% initialize
t0= datetime('now');

% iterate over price schedules
for iter_prices= 1:max_iters_LoM


    % prices: over (s,K,I-)
    R= zeros(ns,nK,nI);
    W= zeros(ns,nK,nI);
    for iI= 1:nI
        for iK= 1:nK
            for is= 1:ns
                % get productivity index
                iz= floor((is-1)/nt)+1;
                % interest rate
                R(is,iK,iI)= alpha* zg(iz) *Kg_agg(iK)^(alpha-1) *Hl0(is,iK,iI)^(1-alpha);
                % wage
                W(is,iK,iI)= (1-alpha) *zg(iz) *Kg_agg(iK)^alpha *Hl0(is,iK,iI)^(-alpha);
            end
        end
    end



    % preallocate
    VV1=    zeros(ni*ns,nk,nK,nI);
    iIpf_4= zeros(ni*ns,nk,nK,nI);
    Ipf_4=  zeros(ni*ns,nk,nK,nI);

    % iterate value function
    for iter= 1:max_iters_vf

        % loop over aggregate capital
        parfor iK= 1:nK
            % note: only for Matlab to stop complaining about fake
            % broadcast variables
            tgg= tg;  VVV0= VV0; PPi= Pi;

            % loop over aggregate investment
            for iI= 1:nI
        
                % policy functions over (i,k)
                Lg= zeros(ni,nk);
                Cg= zeros(ni,nk);

                % loop over exogenous stochastic states
                for is= 1:ns
                    % get productivity & indexes
                    iz= floor((is-1)/nt)+1;
                    it= is - nt*floor((is-1)/nt);
                    % get stacked indexes
                    ii= (is-1)*ni+1:is*ni;
            
                    % labor over (i,k|z,tau)
                    bl=  ( ig - R(is,iK,iI).*kg' )./( (1-tgg(it)) *W(is,iK,iI) );
                    lg=  ( bl + sqrt(bl.^2 + 4) )/2;
                    % consumption over (i,k|z,tau)
                    cg=  ( 1 - tgg(it) ).*W(is,iK,iI).*lg  +  R(is,iK,iI).*kg'  -  ig;
            
                    % stack policy functions
                    Lg(ii,:)= lg;
                    Cg(ii,:)= cg;
            
                end
                    
                % period utility over (i,k)
                U=  log(Cg) - Lg.^2/2;
                
                % preallocate
                V1=     zeros(ni*ns,nk);
                iIpf=   zeros(ni*ns,nk);
                Ipf=    zeros(ni*ns,nk);
             
                
                % loop over exogenous states (z,tau)
                for is= 1:ns
                    % get productivity index
                    iz= floor((is-1)/nt)+1;
    
                    % get indexes for stacking
                    ii= (is-1)*ni+1:is*ni;
                    
                    % take expectation over (z,tau)
                    EV= zeros(ni,nk,nK,nI);
                    for iKp= 1:nK
                        for iIp= 1:nI
                            EVk= EV(ni,nk,iKp,iIp);
                            for isp= 1:ns
                                % compute expected value
                                EVk= EVk + VVV0( (isp-1)*ni+1:isp*ni,:, iKp,iIp )*PPi(is,isp);
                            end
                        % store
                        EV(:,:,iKp,iIp)= EVk;
                        end
                    end
    
                    % set interpolant
                    EVinterpolant= griddedInterpolant(xx,yy,zz,ww,EV,'linear','linear');
                
                    % loop over k
                    for ik= 1:nk
                        
                        % future capital: k'(i-,i;k)
                        kp= (1-delta)*kg(ik) + Phi;
                    
                        % interpolate continuation value over (i-,i)
                        Vcont=  EVinterpolant( ones(ni,1)*ig', kp, ones(ni)*Hk0(is,iK,iI), ones(ni)*Hi0(is,iK,iI) );
                        % set objective function: over (i-,i)
                        Uobj=   (1-beta)*U(ii,ik)' + beta*Vcont;
                        % maximize over i
                        [Vk,iimax]= max( Uobj, [],2 );
                    
                        % store VF & PF: over (i-;k)
                        V1(ii,ik)=   Vk;
                        iIpf(ii,ik)= iimax;
                        Ipf(ii,ik)=  ig(iimax);
                    
                    end
    
                    % store
                    VV1(:,:,iK,iI)=     V1;
                    iIpf_4(:,:,iK,iI)=  iIpf;
                    Ipf_4(:,:,iK,iI)=   Ipf;
                end
            end
        end
            
        % check for convergence
        err= max(max(max(max(abs(VV0-VV1)))));
        if (err<tol_vf), break; end

        % update value function
        VV0= VV1;
        
        % print progress
        t1= datetime('now');
        dt= t1-t0;  
        dt.Format= 'mm:ss.SSS';
        fprintf('iter= %3d  |  err= %.2e  |  etime=%s \n', iter, err,dt );
    end
        

    %----------------------------------
    % compute policy functions
    %----------------------------------
    % rename
    Ipf= Ipf_4;
    iIpf= iIpf_4;

    % preallocate
    Lpf= zeros(ni*ns,nk,nK,nI);
    Cpf= zeros(ni*ns,nk,nK,nI);
    Kpf= zeros(ni*ns,nk,nK,nI);

    % loop over aggregate states
    for iK= 1:nK
        for iI= 1:nI
            % loop over exogenous states
            for is= 1:ns

                % get productivity & tau indexes
                iz= floor((is-1)/nt)+1;
                it= is - nt*floor((is-1)/nt);

                % loop over endogenous states
                for ik= 1:nk
                    for i0= 1:ni
                        % get indexes for stacking
                        j= (is-1)*ni+i0;

                        % labor policy
                        bl= ( Ipf(j,ik,iK,iI) - R(is,iK,iI)*kg(ik) )/( (1-tg(it))*W(is,iK,iI) );
                        Lpf(j,ik,iK,iI)=    ( bl + sqrt(bl.^2 + 4) )/2;
                        % consumption policy
                        Cpf(j,ik,iK,iI)=    (1-tg(it))*W(is,iK,iI)*Lpf(j,ik,iK,iI) + R(is,iK,iI)*kg(ik) - Ipf(j,ik,iK,iI);
                        % future-capital policy
                        Kpf(j,ik,iK,iI)=    (1-delta)*kg(ik) + Phi(i0,iIpf(i0,ik,iK,iI));

                    end
                end
            end
        end
    end


    %----------------------------------
    % laws of motion
    %----------------------------------
    [ii2,kk2]=  ndgrid(ig,kg);
    Hk= zeros(ns,nK,nI);
    Hi= zeros(ns,nK,nI);
    Hl= zeros(ns,nK,nI);
    % loop over aggregate states
    for iK= 1:nK
        for iI= 1:nI
            for is= 1:ns
                % get policy function conditional on (s,K,I)
                ii= (is-1)*ni+1:is*ni;
                Kpf_sKI= Kpf(ii,:,iK,iI);
                Ipf_sKI= Ipf(ii,:,iK,iI);
                Lpf_sKI= Lpf(ii,:,iK,iI);
                
                % interpolate capital policy function
                Hhat=   griddedInterpolant(ii2,kk2,Kpf_sKI);
                Hk(is,iK,iI)=  Hhat( Ig_agg(iI), Kg_agg(iK) );
                % interpolate investment policy function
                Hhat=   griddedInterpolant(ii2,kk2,Ipf_sKI);
                Hi(is,iK,iI)=  Hhat( Ig_agg(iI), Kg_agg(iK) );
                % interpolate labor policy function
                Hhat=   griddedInterpolant(ii2,kk2,Lpf_sKI);
                Hl(is,iK,iI)=  Hhat( Ig_agg(iI), Kg_agg(iK) );
            end
        end
    end

    
    % update laws of motion
    adjfac= 0.15;
    Hk0=    adjfac*Hk + (1-adjfac)*Hk0;
    Hi0=    adjfac*Hi + (1-adjfac)*Hi0;
    Hl0=    adjfac*Hl + (1-adjfac)*Hl0;

    % check convergence
    err_k=  max(max(max(abs(Hk0-Hk))));
    err_i=  max(max(max(abs(Hi0-Hi))));
    err_l=  max(max(max(abs(Hl0-Hl))));
    err=    max([err_k,err_i,err_l]);

    % print progress
    t1= datetime('now');
    dt= t1-t0;  
    dt.Format= 'mm:ss.SSS';
    fprintf('\niter= %3d  |  err(K,I,L)= (%.1e,%.1e,%.1e)  |  etime=%s \n', ...
        iter_prices, err_k,err_i,err_l,dt );
    if (err<tol_LoM), break; end

end



% store results
if (store_results==1)
    save('Hk0.mat','Hk0');
    save('Hi0.mat','Hi0');
    save('Hl0.mat','Hl0');
    save('VV0.mat','VV0');
end





%%
%-------------------------------------------------------------------------------
%  Figures
%-------------------------------------------------------------------------------
figure(1);
tl= tiledlayout('flow');
for is= 1:ns
    ii= (is-1)*ni+1:is*ni;
    iz= floor((is-1)/nt)+1;
    it= is - nt*floor((is-1)/nt);

    % interest rate
    nexttile;
    hold on;
    surf(Kg_agg,Ig_agg,squeeze(R(is,:,:))'); 
    scatter3(kss,iss,rss, 'red', 'filled');
    hold off;
    view(150,30);
    xlabel('capital','rotation',0); ylabel('investment','rotation',0); 
    title(sprintf('$(z=%.2f,\\tau=%.2f)$',zg(iz),tg(it)), 'interpreter','latex');
end
title(tl, 'Interest Rate Schedule');


figure(2);
tl= tiledlayout('flow');
for is= 1:ns
    ii= (is-1)*ni+1:is*ni;
    iz= floor((is-1)/nt)+1;
    it= is - nt*floor((is-1)/nt);

    % investment
    nexttile;
    hold on;
    plot( ig, Ipf( ii,(nk+1)/2, 4,4) );
    plot( ig, ig, 'w--' );
    hold off;
    xline(iss);
    title(sprintf('$(z=%.2f,\\tau=%.2f)$',zg(iz),tg(it)), 'interpreter','latex');
    xlabel('investment');
end
title(tl, '$i(i^{-},k=k^{ss})$', 'interpreter','latex');


figure(3);
tl= tiledlayout('flow');
for is= 1:ns
    ii= (is-1)*ni+1:is*ni;
    iz= floor((is-1)/nt)+1;
    it= is - nt*floor((is-1)/nt);

    % consumption
    nexttile;
    hold on;
    surf(ig,kg, squeeze(Cpf(ii,:,4,4))'); 
    hold off;
    view(-30,30);
    xlabel('investment','rotation',0); ylabel('capital','rotation',0); 
    title(sprintf('$(z=%.2f,\\tau=%.2f)$',zg(iz),tg(it)), 'interpreter','latex');
end
title(tl, 'Consumption');


