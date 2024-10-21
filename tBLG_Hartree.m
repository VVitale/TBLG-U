function [scale_axis1,allbands,all_kpts1,qvecs,vkp,...
    gstar_b1,gstar_b2,drhoG,Ef,tot_dim,Vc] = tBLG_Hartree(n,a1,hvf,knum,...
              gstar_cut_fac0,max_iter,V0,V1,U_param,full_bz,plot_DOS,write_eigvecs,...
              u1, u2, ax_m,valley,dE,bz_n,epsilon,...
              filling,potname,Temp)

max_inter_q = 3; % [#: includes] 3: nearest 
unit_dim=2; % Number of layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining twist angle based on integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(n~=0)
   tar_theta = acos((3*n^2 + 3*n + 0.5)/(3*n^2 + 3*n + 1));
   L = a1*sqrt(3*n^2 + 3*n + 1);
else
   tar_theta = 0;
   L = a1;
end
Vc = L^2*sqrt(3)/2;

tar_theta = tar_theta*180/pi;
fprintf(' Hartree calculation for TBLG with theta=%2.2f° \n', tar_theta);

rot_theta = tar_theta/180*pi;
gstar_b1=2*pi/L*[1/sqrt(3),1];
gstar_b2=4*pi/L*[-1/sqrt(3),0];
gstar_length=norm(gstar_b1);

num_inter_qs=3;

gstar_max=100; % VV: Dimension of wavefunctions exapansion, before applying cutoff

gstar_coor=[];
gstar_cut=gstar_cut_fac0*gstar_length;

%VV: Generate G vectors to include in expansion of wavefunction ?
ind=0;
for ind1=(-gstar_max):gstar_max
    for ind2=(-gstar_max):gstar_max
        vec=gstar_b1*ind1+gstar_b2*ind2;% VV: G vectors
        if norm(vec)<gstar_cut % VV: If length of G is within circle of radius gstar_cut centred at Gamma
            ind=ind+1;
            gstar_coor(ind,1:2)=vec(1:2); % VV: Store coordinates of G vectors
        end   
    end
end

num_gstar=ind; % VV: Number of G vectors contained in circle

fprintf('Number of G vectors: %i\n',num_gstar)

% VV: Dimension of Hamiltonian matrices: # Gvecs * # Layers * # Sublattices
tot_dim=num_gstar*unit_dim*2;

fprintf('Size of Hamiltonian: %i\n',tot_dim)
all_index_L1=reshape(1:(num_gstar*unit_dim),unit_dim,num_gstar);
all_index_L2=reshape((num_gstar*unit_dim+1):(2*num_gstar*unit_dim),unit_dim,num_gstar);
%all_index=reshape(1:tot_dim,unit_dim,tot_dim/unit_dim);

% Generate k-points
[all_kpts1,scale_axis1,knum_tot,vkp,qvecs] = generate_kpoints3(gstar_b1,gstar_b2,...
    1,knum,1,bz_n,plot_DOS);

[connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,...
    tot_num_G12,connect_Mat_pot2] = initialise_workspace(1,gstar_b1,gstar_b2,num_inter_qs,num_gstar,gstar_coor,potname);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Hartree main loop
%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 1;
rho_thr = 10^(-6);
delta = 10;
%delta_V0 = 10;
%delta_V1 = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Reference state
%%%%%%%%%%%%%%%%%%%%%%%%%%
%[eigvals_ref,eigvecs_ref]=build_H_tBLG(u1,u2,1,gstar_b1,gstar_b2,num_gstar,...
%    gstar_coor,all_index_L1,all_index_L2,...
%    rot_theta,hvf,tot_dim,max_inter_q,knum_tot,...
%    all_kpts1,write_eigvecs,epsilon,0,Vc,...
%    connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,...
%    'NONE',0,0);

%CN = find_chem_pot(knum_tot,sort(eigvals_ref),0,Temp,1000,0.00001);

%eigvals_ref = eigvals_ref - CN;
%CN = find_chem_pot(knum_tot,sort(eigvals_ref),0,Temp,1000,0.00001)
%Ef = find_chem_pot(knum_tot,sort(eigvals_ref),filling,Temp,1000,0.00001)

%[ref_1,ref_2,ref_3,ref_4] = compute_rho(eigvals_ref,eigvecs_ref,connect_Mat_hartree,tot_dim,knum_tot,0,CN)
%rhoG_ref = 4*real(sum(ref_1)+sum(ref_2)+sum(ref_3)+sum(ref_4))/6/knum_tot

%[ref0_1,ref0_2,ref0_3,ref0_4] = compute_rho0(eigvals_ref,eigvecs_ref,tot_dim,knum_tot,0,CN);
%drho0_ref = 4*real(ref0_1+ref0_2+ref0_3+ref0_4)/knum_tot;

%might also need a guess for drhoG ... although it should just occur

%%%%%%%%
%up spin - initial calc
%%%%%%%%
[eigvals_up,eigvecs_up]=build_H_tBLG(u1,u2,1,gstar_b1,gstar_b2,num_gstar,...
    gstar_coor,all_index_L1,all_index_L2,...
    rot_theta,hvf,tot_dim,max_inter_q,knum_tot,...
    all_kpts1,write_eigvecs,epsilon,0,Vc,...
    connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,...
    potname,V0,V1);

%%%%%%%%
%dw spin - initial calc
%%%%%%%%
[eigvals_dw,eigvecs_dw]=build_H_tBLG(u1,u2,1,gstar_b1,gstar_b2,num_gstar,...
    gstar_coor,all_index_L1,all_index_L2,...
    rot_theta,hvf,tot_dim,max_inter_q,knum_tot,...
    all_kpts1,write_eigvecs,epsilon,0,Vc,...
    connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,...
    potname,-V0,-V1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Reference state and deltas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbands = [eigvals_up;eigvals_dw];
allbands = sort(allbands,1);
       
CN = find_chem_pot(knum_tot,allbands,0,Temp,1000,0.00001);
allbands = allbands - CN;
eigvals_up = eigvals_up - CN;
eigvals_dw = eigvals_dw - CN;
CN = find_chem_pot(knum_tot,allbands,0,Temp,1000,0.00001);
Ef = find_chem_pot(knum_tot,allbands,filling,Temp,1000,0.00001);

[ref_up_1,ref_up_2,ref_up_3,ref_up_4] = compute_rho(eigvals_up,eigvecs_up,connect_Mat_hartree,tot_dim,knum_tot,0,0);
[ref_dw_1,ref_dw_2,ref_dw_3,ref_dw_4] = compute_rho(eigvals_dw,eigvecs_dw,connect_Mat_hartree,tot_dim,knum_tot,0,0);
rhoG_ref = 2*real(sum(ref_up_1)+sum(ref_up_2)+sum(ref_up_3)+sum(ref_up_4)...
    +sum(ref_dw_1)+sum(ref_dw_2)+sum(ref_dw_3)+sum(ref_dw_4))/6/knum_tot;

[ref_up0_1,ref_up0_2,ref_up0_3,ref_up0_4] = compute_rho0(eigvals_up-CN,eigvecs_up,tot_dim,knum_tot,0,0);
[ref_dw0_1,ref_dw0_2,ref_dw0_3,ref_dw0_4] = compute_rho0(eigvals_dw-CN,eigvecs_dw,tot_dim,knum_tot,0,0);
rho0_ref = 2*real(ref_up0_1+ref_up0_2+ref_up0_3+ref_up0_4+ref_dw0_1+ref_dw0_2+ref_dw0_3+ref_dw0_4)/knum_tot;

[up_1,up_2,up_3,up_4] = compute_rho(eigvals_up,eigvecs_up,connect_Mat_hartree,tot_dim,knum_tot,0,Ef)
[dw_1,dw_2,dw_3,dw_4] = compute_rho(eigvals_dw,eigvecs_dw,connect_Mat_hartree,tot_dim,knum_tot,0,Ef)
rhoG = 2*real(sum(up_1)+sum(up_2)+sum(up_3)+sum(up_4)+sum(dw_1)+sum(dw_2)+sum(dw_3)+sum(dw_4))/6/knum_tot;

[up0_1,up0_2,up0_3,up0_4] = compute_rho0(eigvals_up-CN,eigvecs_up,tot_dim,knum_tot,0,Ef);
[dw0_1,dw0_2,dw0_3,dw0_4] = compute_rho0(eigvals_dw-CN,eigvecs_dw,tot_dim,knum_tot,0,Ef);
rho0 = 2*real(up0_1+up0_2+up0_3+up0_4+dw0_1+dw0_2+dw0_3+dw0_4)/knum_tot;

drhoG = rhoG - rhoG_ref;
drho0 = rho0 - rho0_ref;

%
%this is the magnetic calculations which need to be done ... 
%

[V0p,V1p] = compute_magnetic_potentials(U_param,potname,knum_tot,...
	up0_1,up0_2,up0_3,up0_4,dw0_1,dw0_2,dw0_3,dw0_4,...
        up_1,up_2,up_3,up_4,dw_1,dw_2,dw_3,dw_4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting self-consistent loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%need to define some parameter for the Hubbard interaction ... 
drhoG_new = drhoG;
drho0_new = drho0;
V0p_new = V0p;
V1p_new = V1p;
delta_all = 1;
while iter < max_iter && delta_all > rho_thr% && delta_V0 > rho_thr && delta_V1 > rho_thr
        fprintf('** Iteration %i \n',iter)
        %%%%%%%%
        %spin up
        %%%%%%%%
        [eigvals_up,eigvecs_up] = build_H_tBLG(u1,u2,1,gstar_b1,gstar_b2,num_gstar,...
           gstar_coor,all_index_L1,all_index_L2,rot_theta,...
           hvf,tot_dim,max_inter_q,knum_tot,all_kpts1,write_eigvecs,...
           epsilon,drhoG,Vc,...
           connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,potname,...
           V0p,V1p);
        
        %%%%%%%%%%
        %spin down
        %%%%%%%%%%
        [eigvals_dw,eigvecs_dw] = build_H_tBLG(u1,u2,1,gstar_b1,gstar_b2,num_gstar,...
           gstar_coor,all_index_L1,all_index_L2,rot_theta,...
           hvf,tot_dim,max_inter_q,knum_tot,all_kpts1,write_eigvecs,...
           epsilon,drhoG,Vc,...
           connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,potname,...
           -V0p,-V1p);

        allbands = [eigvals_up;eigvals_dw];
        allbands = sort(allbands,1);
      
        CN = find_chem_pot(knum_tot,allbands,0,Temp,1000,0.00001);
        allbands = allbands - CN;
        eigvals_up = eigvals_up - CN;
        eigvals_dw = eigvals_dw - CN;
        CN = find_chem_pot(knum_tot,allbands,0,Temp,1000,0.00001);
        Ef = find_chem_pot(knum_tot,allbands,filling,Temp,1000,0.00001);

        [ref_1_up,ref_2_up,ref_3_up,ref_4_up] = compute_rho(eigvals_up,eigvecs_up,connect_Mat_hartree,tot_dim,knum_tot,0,0);
        [ref_1_dw,ref_2_dw,ref_3_dw,ref_4_dw] = compute_rho(eigvals_dw,eigvecs_dw,connect_Mat_hartree,tot_dim,knum_tot,0,0);
        rhoG_ref_new = 2*real(sum(ref_1_up)+sum(ref_2_up)+sum(ref_3_up)+sum(ref_4_up)+sum(ref_1_dw)...
            +sum(ref_2_dw)+sum(ref_3_dw)+sum(ref_4_dw))/6/knum_tot
        
        [up_1,up_2,up_3,up_4] = compute_rho(eigvals_up,eigvecs_up,connect_Mat_hartree,tot_dim,knum_tot,0,Ef)
        [dw_1,dw_2,dw_3,dw_4] = compute_rho(eigvals_dw,eigvecs_dw,connect_Mat_hartree,tot_dim,knum_tot,0,Ef)
        rhoG_new = 2*real(sum(up_1)+sum(up_2)+sum(up_3)+sum(up_4)+sum(dw_1)+sum(dw_2)+sum(dw_3)+sum(dw_4))/6/knum_tot
        
        drhoG_new = rhoG_new - rhoG_ref_new
        
        [ref_up0_1,ref_up0_2,ref_up0_3,ref_up0_4] = compute_rho0(eigvals_up,eigvecs_up,tot_dim,knum_tot,0,0);
        [ref_dw0_1,ref_dw0_2,ref_dw0_3,ref_dw0_4] = compute_rho0(eigvals_dw,eigvecs_dw,tot_dim,knum_tot,0,0);
        rho0_ref_new = 2*real(ref_up0_1+ref_up0_2+ref_up0_3+ref_up0_4+ref_dw0_1+ref_dw0_2+ref_dw0_3+ref_dw0_4)/knum_tot;

        [up0_1,up0_2,up0_3,up0_4] = compute_rho0(eigvals_up,eigvecs_up,tot_dim,knum_tot,0,Ef);
        [dw0_1,dw0_2,dw0_3,dw0_4] = compute_rho0(eigvals_dw,eigvecs_dw,tot_dim,knum_tot,0,Ef);
        rho0_new = 2*real(up0_1+up0_2+up0_3+up0_4+dw0_1+dw0_2+dw0_3+dw0_4)/knum_tot;

        drho0_new = rho0_new - rho0_ref_new
        
        [V0p_new,V1p_new] = compute_magnetic_potentials(U_param,potname,knum_tot,...
		up0_1,up0_2,up0_3,up0_4,dw0_1,dw0_2,dw0_3,dw0_4,...
                up_1,up_2,up_3,up_4,dw_1,dw_2,dw_3,dw_4);

        delta = abs(drhoG_new - drhoG);
        delta_V0 = abs(V0p_new - V0p);
        delta_V1 = abs(V1p_new - V1p);
        delta_all = delta + delta_V0 + delta_V1;
        fprintf(' Delta: %2.8f  Delta_0: %2.8f  Delta_1: %2.8f \n',delta,delta_V0,delta_V1)
        %fprintf(' Im/Re: %2.8f \n',ratio)
        iter = iter + 1;
        
        %how come there is no mixing going on here? 
        mix = 0.1;
        drhoG = (1-mix)*drhoG + mix*drhoG_new;
        drho0 = (1-mix)*drho0 + mix*drho0_new;
        V0p = (1-mix)*V0p + mix*V0p_new;
        V1p = (1-mix)*V1p +mix*V1p_new
end

%Plot bands as surfaces over the entire miniBZ
if (full_bz == 1)
    band_mesh_h = zeros(bz_n,bz_n);
    band_mesh_e = zeros(bz_n,bz_n);

    kx_mesh = zeros(bz_n,bz_n);
    ky_mesh = zeros(bz_n,bz_n);

    nb = size(allbands,1);
    bz_n = sqrt(size(allbands,2));
    idx = 1;
    for y = 1:bz_n
        for x = 1:bz_n
            band_mesh_h(x,y) = allbands(nb/2,idx);
            band_mesh_e(x,y) = allbands(nb/2+1,idx);
            kx_mesh(x,y) = all_kpts1(idx,1);
            ky_mesh(x,y) = all_kpts1(idx,2);
            idx = idx+1;
        end
    end

    clf
    hold on
    surf(kx_mesh,ky_mesh,band_mesh_h,'EdgeColor','none')
    surf(kx_mesh,ky_mesh,band_mesh_e,'EdgeColor','none')
end

if(plot_DOS)
   plot_DOS_fun(allbands,dE,knum_tot,Ef)
end
    
    
figure
for chi = valley
    %
    %can we also do spin?
    %
    % Re-generate k-points for plotting bandstructure
    [all_kpts_BS,scale_axis_BS,knum_tot_BS,~,~] = generate_kpoints3(gstar_b1,gstar_b2,...
        chi,knum,0,0,0);
    
    [connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,tot_num_G12,connect_Mat_pot2] = initialise_workspace(chi,...
       gstar_b1,gstar_b2,num_inter_qs,num_gstar,gstar_coor,potname);
    
   [allbands_BS,~] = build_H_tBLG(u1,u2,chi,gstar_b1,gstar_b2,num_gstar,...
       gstar_coor,all_index_L1,all_index_L2,rot_theta,...
       hvf,tot_dim,max_inter_q,knum_tot_BS,all_kpts_BS,...
       write_eigvecs,epsilon,drhoG,Vc,...
       connect_Mat_L12,connect_Mat_hartree,connect_Mat_pot,potname,...
       V0p,V1p);
   
    if(chi==1)
        plot(scale_axis_BS,allbands_BS-CN,'b','LineWidth',2.0)
        axis([0 1 -ax_m ax_m])
        ylabel('Energy (eV)')
        title(['$\theta = ' num2str(tar_theta,'%.2f') '^\circ$'])
        grid 'on'
        set(gca,'XTick',[0,scale_axis_BS(knum),scale_axis_BS(2*knum),1]);
        set(gca,'XTickLabel',{'K','$\Gamma$','M',"K'"});
        hold on
    elseif(chi==-1)
        plot(scale_axis_BS,allbands_BS-CN,'r','LineWidth',2.0)
        hold on
    end
    box 'on'
    outfname = join(['bands_theta_',num2str(tar_theta,'%.2f'),'_valley_',num2str(chi),'_nu=',num2str(filling),'.dat']);
    fileID = fopen(outfname,'w');
    for in = 1 : tot_dim
      for ik = 1 : knum_tot_BS
         fprintf(fileID,'%2.8f  %2.8f\n',scale_axis_BS(ik),allbands_BS(in,ik)-CN);
      end
         fprintf(fileID,'\n');
    end
    fclose(fileID);
end

filename = outfname;%'kp_band_Hartree_data';

save(filename, ...
    'tar_theta','allbands_BS','scale_axis_BS', ...
    'all_kpts_BS');

clear all_kpts_BS scale_axis_BS knum_tot_BS allbands_BS ...
      gstar_coor vec all_index_L1 all_index_L2 all_index

end
