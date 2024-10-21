function [allbands1,eigvecs]= build_H_tBLG(u1,u2,chi,gstar_b1,gstar_b2,num_gstar,...
    gstar_all,all_index_L1,all_index_L2,...
    rot_theta,Dirac_v1,tot_dim,max_inter_q,...
    knum_tot,all_kpts1,write_eigvecs,epsilon,drhoG,Vc,...
    cmat_interlayer,cmat_hartree,cmat_pot,potential,...
    V0,V1)
    %drhoG
    %V0
    %V1

    % Units factors
    ang2bohr = 1.8897259886;
    hartree2eV = 27.211396641308;

    % Generate U matrices in K-M PRB
    interlayer_pot=zeros(2,2,3);
    
    interlayer_pot(:,:,1)=[u1,u2;u2,u1];
    interlayer_pot(:,:,2)=[u1,u2*exp(-1i*chi*2*pi/3);u2*exp(1i*chi*2*pi/3),u1];
    interlayer_pot(:,:,3)=[u1,u2*exp(1i*chi*2*pi/3);u2*exp(-1i*chi*2*pi/3),u1];

    % Generate Hartree contribution diagonal matrix 
    % Metallic gates equally distant from bilayer
    L = sqrt(2/sqrt(3)*Vc)*ang2bohr;
    hartree_pot = 1/L/epsilon*drhoG*hartree2eV*eye(2);
    if(~strcmp(potential,'NONE'))
       if(strcmp(potential,'HAFM'))
          mag_pot = 1i*V1*[1 0;0 -1];
       elseif(strcmp(potential,'NAFM'))
          mag_pot = V1*[1 0;0 -1];
       elseif(strcmp(potential,'MAFM'))
          mag_pot = V0*[1 0;0 -1];%eye(2);% LOOKS WRONG TO ME
          mag_pot1 = V1*[1 0;0 -1];
       elseif(strcmp(potential,'FM'))
          mag_pot = V1*eye(2);%might need a spin index here ... 
          %mag_pot1 = V1*eye(2);
          %mag_pot2 = 1i*V2*eye(2);
       end
    else
       mag_pot = zeros(2);
    end
    %*tanh(norm(gstar_b1)*gate_dis)

    % Spin rotations
    expfac1=exp(-1i*chi*rot_theta/2);
    expfac2=exp(1i*chi*rot_theta/2);
    
    % K and K'
    K1 = chi*(-gstar_b1-2*gstar_b2)/3;
    K2 = chi*(gstar_b1-gstar_b2)/3;

    % Rotated Dirac Hamiltonian function
    Dirac_ham_L1=@(kknow) -Dirac_v1*[0,expfac1*(1i*kknow(2)-chi*kknow(1)+chi*K1(1)-1i*K1(2));conj(expfac1)*(-1i*kknow(2)-chi*kknow(1)+chi*K1(1)+1i*K1(2)),0];
    Dirac_ham_L2=@(kknow) -Dirac_v1*[0.0,expfac2*(1i*kknow(2)-chi*kknow(1)+chi*K2(1)-1i*K2(2));conj(expfac2)*(-1i*kknow(2)-chi*kknow(1)+chi*K2(1)+1i*K2(2)),0.0];
    
    Hmat_inter=zeros(tot_dim);

    % Interlayer TBLG Hamiltonian only bottom triangular part
    for indq=1:max_inter_q
        T_tmp=squeeze(interlayer_pot(:,:,indq)); % VV: U matrices in K-M PRB paper
        Hmat_inter(all_index_L2(:),all_index_L1(:)) = Hmat_inter(all_index_L2(:),all_index_L1(:))+kron(squeeze(cmat_interlayer(:,:,indq)),T_tmp);
        % VV: Here he's using the kroneker tensor product to generate all off diagonal Hamiltonian matrix
    end
    
    % Add upper triangular part as well
    Hmat_inter=Hmat_inter+Hmat_inter';

    allbands1=zeros(tot_dim,knum_tot);
    eigvecs=zeros(tot_dim,tot_dim,knum_tot);
    
    % Add Hartree contribution to each layer
    Hmat_hartree = zeros(tot_dim);
    for indq=1:2*max_inter_q
        Hmat_hartree(all_index_L1(:),all_index_L1(:)) = Hmat_hartree(all_index_L1(:),all_index_L1(:))+kron(squeeze(cmat_hartree(:,:,indq)),hartree_pot);
        Hmat_hartree(all_index_L2(:),all_index_L2(:)) = Hmat_hartree(all_index_L2(:),all_index_L2(:))+kron(squeeze(cmat_hartree(:,:,indq)),hartree_pot);
    end
    % Add potential contribution to each layer
    if(~strcmp(potential,'NONE'))
       Hmat_pot = zeros(tot_dim);
       if(strcmp(potential,'HAFM'))
          for indq=1:2*max_inter_q
              Hmat_pot(all_index_L1(:),all_index_L1(:)) = Hmat_pot(all_index_L1(:),all_index_L1(:))+kron(squeeze(cmat_pot(:,:,indq)),chi*mag_pot);
              Hmat_pot(all_index_L2(:),all_index_L2(:)) = Hmat_pot(all_index_L2(:),all_index_L2(:))+kron(squeeze(cmat_pot(:,:,indq)),-chi*mag_pot);
          end
       elseif(strcmp(potential,'NAFM'))
          for indq=1:2*max_inter_q
              Hmat_pot(all_index_L1(:),all_index_L1(:)) = Hmat_pot(all_index_L1(:),all_index_L1(:))+kron(squeeze(cmat_pot(:,:,indq)),mag_pot);
              Hmat_pot(all_index_L2(:),all_index_L2(:)) = Hmat_pot(all_index_L2(:),all_index_L2(:))+kron(squeeze(cmat_pot(:,:,indq)),mag_pot);
          end
       elseif(strcmp(potential,'MAFM'))
          Hmat_pot(all_index_L1(:),all_index_L1(:)) = Hmat_pot(all_index_L1(:),all_index_L1(:))+kron(squeeze(cmat_pot(:,:,1)),mag_pot);
          Hmat_pot(all_index_L2(:),all_index_L2(:)) = Hmat_pot(all_index_L2(:),all_index_L2(:))+kron(squeeze(cmat_pot(:,:,1)),mag_pot);
          for indq=2:2*max_inter_q + 1%why does this one start from 2 and go to +1?
              Hmat_pot(all_index_L1(:),all_index_L1(:)) = Hmat_pot(all_index_L1(:),all_index_L1(:))+kron(squeeze(cmat_pot(:,:,indq)),mag_pot1);
              Hmat_pot(all_index_L2(:),all_index_L2(:)) = Hmat_pot(all_index_L2(:),all_index_L2(:))+kron(squeeze(cmat_pot(:,:,indq)),mag_pot1);
          end
         %
         %there is nothing for the FM order?!
         %
       end
    end

    % Start diagonalisation
    for indk=1:knum_tot
        %fprintf("%d/%d k sampling \n",indk,knum_tot);
        know=all_kpts1(indk,1:2);
        shift_klist_L1=zeros(num_gstar,2);
        shift_klist_L2=zeros(num_gstar,2);
        shift_klist_L1(:,1)=gstar_all(:,1)+know(1);
        shift_klist_L1(:,2)=gstar_all(:,2)+know(2);
        shift_klist_L2(:,1)=-gstar_all(:,1)+know(1);
        shift_klist_L2(:,2)=-gstar_all(:,2)+know(2);
        
        Hmat=Hmat_inter+Hmat_hartree;
        if(~strcmp(potential,'NONE'))
           Hmat = Hmat + Hmat_pot;
        end
        % intralayer part
        for indh=1:num_gstar
            Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Dirac_ham_L1(shift_klist_L1(indh,:));
            Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Dirac_ham_L2(shift_klist_L2(indh,:));
        end
        [V,D] = eig(Hmat,'vector');
        [allbands1(:,indk),indband]=sort(real(D),'ascend');
        eigvecs(:,:,indk) = V(:,indband);
        if(write_eigvecs)
            for indh=1:num_gstar
                if(indh<=size(allbands1,1)/2)
                    vec_outfname=join(['val_eigvecs.',num2str(indk),'.',num2str(indh),'.dat']);
                    file_ID=fopen(vec_outfname,'w');
                    eig_outfname=join(['val_eigvals.',num2str(indk),'.dat']);
                    file_ID2=fopen(eig_outfname,'w');
                else
                    vec_outfname=join(['cond_eigvecs.',num2str(indk),'.',num2str(indh),'.dat']);
                    file_ID=fopen(vec_outfname,'w');
                    eig_outfname=join(['cond_eigvals.',num2str(indk),'.dat']);
                    file_ID2=fopen(eig_outfname,'w');
                end
                fprintf(file_ID,"%s\n",'# First line is a comment');
                fprintf(file_ID2,"%s\n",'# First line is a comment');
                fprintf(file_ID," %15.8f     %15.8f\n",real(eigvecs(:,indh,indk)),imag(eigvecs(:,indh,indk)));
                fprintf(file_ID2," %15.8f\n",allbands1(:,indk));
                fclose(file_ID);
                fclose(file_ID2);
            end
        end
    end
clear D V Hmat Hmat_inter Hmat_hartree interlayer_pot hartree_pot interall_given_qs ...
      hartree_given_qs interlayer_qvecs hartree_qvecs cmat_interlayer
end
