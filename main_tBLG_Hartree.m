%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                            %%
%% Continuum model for twisted bilayer graphene with atomic relaxation.       %%
%% This routine computes the band structure of twisted bilayer graphene with  %%
%% self-consistent Hartree theory.                                            %%
%% Also incorporates local Hubbard interactions                               %%
%%                                                                            %%
%% Ref.                                                                       %%
%% [1] R. Bistritzer,and A. MacDonald, PNAS 108 (30) 12233-12237 (2011)       %%
%% [2] T. Cea, N. R. Walet, and F. Guinea, PRB 100, 205113 (2019)             %%
%% [3] A. Jimeno-Pozo, Z.A.H. Goodwin, V. Vitale, et al. Adv. Phys. Res., 2:  %% 
%% 2300048 (2023)                                                             %%
%% [4] https://github.com/stcarr/kp_tblg                                      %%
%% [5] S. Fang, S. Carr et al. arXiv:1908.00058 (2019)                        %%
%% [6] S. Carr, S. Fang et al. arXiv:1901.03420 (2019)                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                            %%
%% Written by Valerio Vitale                                                  %%
%% - modified by ZAHG                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT VARIABLES
% Integer for twist angle. Here we use the convention m=n+1
n = 31;
% Lattice constant
a1 = 2.46;
% Fermi velocity of graphene in eV*Ang
hbar_vf = 2.1354*a1;
% Number of kpoints along the first segment of the k-path
knum = 15;
% Cut-off for G-vector expansion
cut_fac0 = 5;
% Value of u1 and u2 as in Ref. [2]
u1 = 0.0797;
u2 = 0.0975;
% Valley: +1 for K and -1 for Kprime
valley = [-1,1];
% Max number of iterations in SCF Hartree
max_iter = 5;
% Filling of flat bands. 0:CNP, -4:all flat valence bands are empty, 
% +4:all conduction bands are filled. If filling is different from 0
% a Hartree calculation is performed
filling = 0;
% Value of the dielectric constant
epsilon = 20;
% Whether to add a magnetic potential and which one
potential = 'HAFM';

%these are now the guess values to prepare the system with that ordering
V0 = 500e-5;
V1 = 100e-5;
%V2 = 5e-6;

% values for plotting
dE = 0; % For plotting DOS in eVvi
bz_n = 12; % Square root of total number of k-points in uniform grid 
full_bz = 0;
plot_DOS = 0;
% Whether to write eigenvectors as text file
write_eigvecs = false;
ax_m = 0.05;% Delta energy around Fermi level for plotting in eV 

% meV Hubbard parameter 
U_param = 0.01;%10
% Electronic temperature in K for Fermi level
e_temp = 0.1;

%
%can set up some loops to go over different fillings ... 
%or different twist angles ... 
%

[scale_axis1,allbands1,all_kpts1,qvecs,vkp,G1,G2,drhoG,Ef,tot_dim,Vc] = tBLG_Hartree(n,a1,hbar_vf,knum,...
              cut_fac0,max_iter,V0,V1,U_param,full_bz,plot_DOS,write_eigvecs,...
              u1, u2, ax_m,valley,dE,bz_n,epsilon,filling,potential,e_temp);

%
%then write some stuff to print eigenvalues ... or total energies ... 
%
