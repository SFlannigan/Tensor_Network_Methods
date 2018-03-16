%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MPS Ground State Calculation from Product State
%%  Using MPO-TDVP algorithm.
%%  - With particle number conservation
%%
%%
%%  N Identical Bosons in an M Site Lattice
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all
addpath('../Kernel/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters
M=6; % Number of lattice sites
N =3;% Total number of particles
N_max =3; % Maximum number of particles allowed per site
Bond_Dimension=10; % Maximum mps Bond Dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state=mps_var(M,Bond_Dimension,N,N_max);
%%%%%%%%%%%%%%%%%%%
% Define Initial Particle Distribution. Vector length must be the same as
% the number of particles. If this step is not done, the algorithm uses a
% random distribution instead.
state=state.set_Particle_Position([2,2,4]); 
%%%%%%%%%%%%%%%%%%%
state=state.set_rand_product_state_for_Variational_Algorithms(Bond_Dimension);
state = state.Canonicalise_1s('R-L','false');

% Check Normalisation
Check_Norm_In = state.Full_Norm;

% Find initial particle distribution
Num_2 = zeros(1,M);
for site = 1:M
    Num_2(site) = state.Site_Site_Particle_Corr(site,site);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TDVP Ground State Calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate MPO Hamiltonian 
J=1; U=10; E=0*ones(M,1); u_chem=0;
H=mpo_cpn(M,N_max);
H=H.Simple_1D_Nearest(J,U,E-u_chem);
% J=-1;Jdash=-sqrt(2);
% H=H.SawTooth_1D(J,Jdash,U,E-u_chem);

% Increase Maximum bond dimension. Adds Zeros to the MPS entries.
% Bond_Dim = 50;
% state=state.Increase_bond_dim(Bond_Dim); 

% Time Evolution Parameters
T=1;
dt=0.1;

state=state.TDVP_2s(H,dt,T,'false');

% Check Normalisation
Check_Norm = state.Full_Norm;

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

