%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MPS Ground State Calculation from Product State
%%  Using Imaginary time evolution algorithm (TEBD)
%%  - With particle number conservation
%%
%%  Bosons in Y-Junction
%%  
%%  Two Right Branches are connected by a long range interaction term 
%%  implemented with SWAP gates
%%
%%  MPS Layout:
%%  LEFT_BRANCH - CENTRAL_SITE - RIGHT_BRANCH_1 - RIGHT_BRANCH_2
%%
%% Many SWAP gates are implemented in the same operation to connect the 
%% Central site with the second right branch.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all
addpath('../Kernel/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters
M_L=5; % Number of lattice sites to left of junction
M_R1 = 5; % Number of lattice sites in arm one to right of junction
M_R2 = 5; % Number of lattice sites in arm two to right of junction
M=M_L+M_R1+M_R2+1;
N =3;% Total number of particles
N_max =2; % Maximum number of particles allowed per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state=mps_cpn(M,1,N,N_max);
% state=state.set_Particle_Position([12,16,16]);
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
%% TEBD (imaginary) time-evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hamiltonian Parameters
J=1; U=1; E=0*ones(M,1); u_chem=0;

x = linspace(1,M,M)'; x_centre=floor(M_L/2); x_un = 1;
E = -2*exp(-(x-x_centre).^2/(2*x_un^2));

% Time-evolution parameters
dt = -0.05i;
T=20;
time_steps=abs(T/dt);

% Choose Maximum bond dimension. The MPS in the future will be truncated to
% this value. This routine will not immediately change the size of the
% tensors, but allow the tensors to naturally grow to this value.
truncation = 5;
state=state.set_bond_dim(truncation); 

% Set vector containing Suzuki-Trotter time steps for TEBD time evolution.
order = 4;
state = state.set_Suzuki_Trotter_order(order);

Error=0;
for tt = 1:time_steps
    tic
    
    [state,Total_error]=state.TEBD_Local_2s_Gates_Y_Junction(dt,J,U,E-u_chem,M_L,M_R1,M_R2);
    Error = Error + Total_error;
    
    % Check Normalisation
    Check_Norm = state.Full_Norm;
    
    % Particle Number
    Num=0;
    for site = 1:M
        Num = Num+state.Site_Site_Particle_Corr(site,site);
    end
    [state,Total_error] = state.Canonicalisation_2s('L-R');
    Energy = state.Get_Energy_2s_Y_Junction(J,U,E,M_L,M_R1,M_R2);
    Var = state.Get_Var_2s_Y_Junction(J,U,E,M_L,M_R1,M_R2);
    
    disp(['step #' num2str(tt)...
        ' -- Variance=' num2str(Var)...
        ' -- Energy=' num2str(Energy)...
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
        Corr(site1,site2) = abs(state.Site_Site_Particle_Corr(site1,site2));
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


