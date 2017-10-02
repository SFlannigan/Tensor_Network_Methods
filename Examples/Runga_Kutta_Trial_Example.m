%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MPS Ground State Calculation from Product State
%%  Using Runga-Kutta like algorithm for time evolution.
%%  - With particle number conservation
%%
%%  [1] - Juan Jos» GarcÃa-Ripoll, 'Time evolution of Matrix Product States'
%%  IOP, (2006)
%%
%%  N Identical Bosons in an M Site Lattice
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all
addpath('../Kernel/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters
M=10; % Number of lattice sites
N =2;% Total number of particles
N_max =N; % Maximum number of particles allowed per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state=mps_cpn(M,1,N,N_max);
%%%%%%%%%%%%%%%%%%%
% Define Initial Particle Distribution. Vector length must be the same as
% the number of particles. If this step is not done, the algorithm uses a
% random distribution instead.
state=state.set_Particle_Position([2,3]); 
%%%%%%%%%%%%%%%%%%%
state=state.set_rand_product_state;
state=state.set_bond_dim(1); % Change bond dimension
[state,Total_error] = state.Canonicalisation_2s('L-R');

% Check Normalisation
Check_Norm = state.Full_Norm;

% Find initial particle distribution
Num_2 = zeros(1,M);
for site = 1:M
    Num_2(site) = state.Site_Site_Particle_Corr(site,site);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runga-Kutta (imaginary) time-evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate MPO Hamiltonian 
J=1; U=1; E=0*ones(M,1); u_chem=0;
H=mpo_cpn(M,N_max);
H=H.Simple_1D_Nearest(J,U,E-u_chem);

% Time-evolution parameters
dt = -0.05i;
T=2;
time_steps=abs(T/dt);

% Choose Maximum bond dimension. The MPS in the future will be truncated to
% this value. This routine will not immediately change the size of the
% tensors, but allow the tensors to naturally grow to this value.
truncation = 10;
state=state.set_bond_dim(truncation); 

% Set vector containing roots of exp(-i*dt*H).
order = 4;
state = state.Set_Runga_Kutta_Order(order);

Error=0;
for tt = 1:time_steps
    tic
    
    [state,Total_error]=state.Runga_Kutta_Time_Evolve(H,dt);
    Error = Error + Total_error;
    
    % Check Normalisation
    Check_Norm = state.Full_Norm;
    
    % Particle Number
    Num=0;
    for site = 1:M
        Num = Num+state.Site_Site_Particle_Corr(site,site);
    end
    
    disp(['step #' num2str(tt)...
        ' -- time=' num2str((tt)*dt)...
        ' -- cpu time=' num2str(toc)...
        ': P_Num=' num2str(Num)...
        ': Truncation Error=',num2str(Error)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Observables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find particle site-site correlation function
Corr = zeros(M,M);
for site1 = 1:M
    for site2 = 1:M
        Corr(site1,site2) = state.Site_Site_Particle_Corr(site1,site2);
    end
end
state.plot_Corr(Corr);
title(['TEBD Site-Site Correlation Function: J=',num2str(J),'; U=',num2str(U),'       '],'fontsize',20);

% Find particle site population function
Num = zeros(1,M);
for site = 1:M
    Num(site) = state.Site_Site_Particle_Corr(site,site);
end
state.plot_Corr(Num);
title(['TEBD Site Population: J=',num2str(J),'; U=',num2str(U),'       '],'fontsize',20);


