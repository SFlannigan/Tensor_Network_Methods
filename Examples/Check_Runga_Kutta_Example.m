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
addpath('../Kernel/E_D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters
M=4; % Number of lattice sites
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

B = Basis_set(N,M);
W_In=state.Calc_State_Vector(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runga-Kutta (imaginary) time-evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Trial MPO Hamiltonian 
H=mpo_cpn(M,N_max);

% Trial Hopping between site and site+1
site=2;
H=H.Trial(site);

% Time-evolution parameters
dt = 1;
T=1;
time_steps=abs(T/dt);

% Choose Maximum bond dimension. The MPS in the future will be truncated to
% this value. This routine will not immediately change the size of the
% tensors, but allow the tensors to naturally grow to this value.
truncation = 10;
state=state.set_bond_dim(truncation); 

state = state.Apply_MPO(H);
                
% Truncate and Canonicalise
% [state,error] = Canonicalisation_2s(state,'L-R');


W_out=state.Calc_State_Vector(B);

% Check Normalisation
Check_Norm = state.Full_Norm;

% Particle Number
Num=zeros(1,M);
for site = 1:M
    Num(site) = state.Site_Site_Particle_Corr(site,site);
end


