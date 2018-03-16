%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for MPS. Run after 'DMRG_Trial_Example'
%%
%% Compares the wavefunction of the TEBD to Exact Diagonalisation
%%
%%  Don't use for large M
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M>11
    error('ERROR: system is too large for ED')
elseif isempty(state.P_Positions)
%     error('ERROR: Must define initial particle distribution with function "mps_cpn.set_Particle_Position".')
end

addpath('../Kernel/E_D');

B = Basis_set(N,M);
H = Onsite_Ham(B,E)+Hop_Ham_SawTooth(B,J,Jdash)+Int_Ham(B,U*ones(1,M),U*ones(1,M));

P_Dist = state.P_Positions;
B_Check = zeros(1,size(B,2));

for ii = 1:size(P_Dist,2)
    B_Check(P_Dist(ii))=B_Check(P_Dist(ii))+1;
end

for ii = 1:size(B,1)
    if all(B(ii,:)==B_Check)
        ind = ii;
    end
end

[Psi,Energy]=eigs(H,1,'smallestreal');

% Psi=Psi/sqrt(Psi'*Psi);
[Num_Check] = Get_Site_Occupation(Psi,B);
W_TEBD = state.Calc_State_Vector(B);
% W_TEBD=W_TEBD/sqrt(W_TEBD'*W_TEBD);

Error = 1 - abs(W_TEBD'*Psi)^2
% Num'-Num_Check
