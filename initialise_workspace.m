function [cmat_inter,cmat_hartree,cmat_pot,tot_num_Ginter,cmat_pot2] = initialise_workspace(chi,gstar_b1,gstar_b2,num_inter_gstar,num_gstar,gstar_all,potname)

    potname
    % Define theshold
    vecthres=1E-5;
    tot_num_Ginter=num_inter_gstar;
    

    % These are all of the Gvecs for transitions in TBLG Hamiltonian
    interlayer_gstar_vecs(1,:)=[0,0];
    interlayer_gstar_vecs(2,:)=chi*gstar_b1;
    interlayer_gstar_vecs(3,:)=chi*(gstar_b1+gstar_b2);%+gstar_b2);
    % Creates matrices out of all those transitions
    inter_qvecs=interlayer_gstar_vecs(:,1:2);
    cmat_inter=zeros(num_gstar,num_gstar,tot_num_Ginter);

    % These are all of the Gvecs for transitions in Hartree Hamiltonian
    hartree_gstar_vecs(1,:)=chi*gstar_b1;
    hartree_gstar_vecs(2,:)=chi*gstar_b2;
    hartree_gstar_vecs(3,:)=chi*(gstar_b1+gstar_b2);
    % Creates matrices out of all those transitions
    hartree_qvecs(1:3,:)=hartree_gstar_vecs(:,1:2);
    hartree_qvecs(4:6,:)=-hartree_gstar_vecs(:,1:2);
    cmat_hartree=zeros(num_gstar,num_gstar,2*tot_num_Ginter);

    if(strcmp(potname,'HAFM') || strcmp(potname,'NAFM'))
       mag_pot_gstar_vecs(1,:)=chi*gstar_b1;
       mag_pot_gstar_vecs(2,:)=chi*gstar_b2;
       mag_pot_gstar_vecs(3,:)=-chi*(gstar_b1+gstar_b2);
       % Creates matrices out of all those transitions
       mag_pot_qvecs(1:3,:)=mag_pot_gstar_vecs(:,1:2);
       mag_pot_qvecs(4:6,:)=-mag_pot_gstar_vecs(:,1:2);
       cmat_pot=zeros(num_gstar,num_gstar,2*tot_num_Ginter);
    elseif(strcmp(potname,'MAFM')) 
       mag_pot_gstar_vecs(1,:)=chi*[0,0];
       mag_pot_gstar_vecs(2,:)=chi*gstar_b1;
       mag_pot_gstar_vecs(3,:)=chi*gstar_b2;
       mag_pot_gstar_vecs(4,:)=chi*(gstar_b1+gstar_b2);
       % Creates matrices out of all those transitions
       mag_pot_qvecs(1,:)=mag_pot_gstar_vecs(1,1:2);
       mag_pot_qvecs(2:4,:)=mag_pot_gstar_vecs(2:4,1:2);
       mag_pot_qvecs(5:7,:)=-mag_pot_gstar_vecs(2:4,1:2);
       cmat_pot=zeros(num_gstar,num_gstar,2*tot_num_Ginter+1);
    elseif(strcmp(potname,'FM'))
       mag_pot_gstar_vecs(1,:)=chi*[0,0];
       mag_pot_gstar_vecs(2,:)=chi*gstar_b1;
       mag_pot_gstar_vecs(3,:)=chi*gstar_b2;
       mag_pot_gstar_vecs(4,:)=chi*(gstar_b1+gstar_b2);
       % Creates matrices out of all those transitions
       mag_pot_qvecs(1,:)=mag_pot_gstar_vecs(1,1:2);
       mag_pot_qvecs(2:4,:)=mag_pot_gstar_vecs(2:4,1:2);
       mag_pot_qvecs(5:7,:)=-mag_pot_gstar_vecs(2:4,1:2);
       cmat_pot=zeros(num_gstar,num_gstar,2*tot_num_Ginter+1);
       cmat_pot2=zeros(num_gstar,num_gstar,2*tot_num_Ginter);
    else
       cmat_pot=0;
       cmat_pot2=0;
    end

    if(~exist('cmat_pot'))
       cmat_pot = 0;
    end
    if(~exist('cmat_pot2'))
       cmat_pot2 = 0;
    end

    % Start loop over G-points and create connecting matrices
    for ind1=1:num_gstar
        tvec1=gstar_all(ind1,:); %VV: G_i
        
        for ind2=1:num_gstar
            tvec2=gstar_all(ind2,:); %VV: -G_i
            
            % VV: Find all q : |q +/- G|  < Gcut
            for indty=1:num_inter_gstar
                tmp_qvec=inter_qvecs(indty,1:2); %VV: q vectors
                dvec_inter=-tvec2-tvec1-tmp_qvec;
                if sqrt(dot(dvec_inter,dvec_inter))<vecthres
                    cmat_inter(ind2,ind1,indty)=1;
                end
            end
            
            for indty=1:2*num_inter_gstar
                tmp_qvec=hartree_qvecs(indty,1:2); %VV: q vectors
                dvec_hartree=-(tvec1)+(tvec2)-tmp_qvec;
                if sqrt(dot(dvec_hartree,dvec_hartree))<vecthres
                   cmat_hartree(ind2,ind1,indty)=1;%2*pi/epsilon/normG*tanh(normG*gate_dis)/Vc;
                end
            end
            if(~strcmp(potname,'NONE'))
               if(strcmp(potname,'NAFM'))
                  for indty=1:2*num_inter_gstar
                      tmp_qvec=mag_pot_qvecs(indty,1:2); %VV: q vectors
                      dvec_pot=-(tvec1)+(tvec2)-tmp_qvec;
                      if sqrt(dot(dvec_pot,dvec_pot))<vecthres
                         cmat_pot(ind2,ind1,indty)=1;
                      end
                  end
               elseif(strcmp(potname,'MAFM') || strcmp(potname,'FM'))
                  for indty=1:2*num_inter_gstar+1
                      tmp_qvec=mag_pot_qvecs(indty,1:2); %VV: q vectors
                      dvec_pot=-(tvec1)+(tvec2)-tmp_qvec;
                      if sqrt(dot(dvec_pot,dvec_pot))<vecthres
                         cmat_pot(ind2,ind1,indty)=1;
                         if(strcmp(potname,'FM'))
                            if(indty <= 3)
                                cmat_pot2(ind2,ind1,indty)=1;
                            else                 
                                cmat_pot2(ind2,ind1,indty)=-1;
                            end
                         end
                      end
                  end
               elseif(strcmp(potname,'HAFM'))
                  for indty=1:2*num_inter_gstar
                      tmp_qvec=mag_pot_qvecs(indty,1:2); %VV: q vectors
                      dvec_pot=-(tvec1)+(tvec2)-tmp_qvec;
                      if sqrt(dot(dvec_pot,dvec_pot))<vecthres
                        if(indty <= 3)
                            cmat_pot(ind2,ind1,indty)=1;
                        else
                            cmat_pot(ind2,ind1,indty)=-1;
                        end
                      end
                  end
               end
            end
        end
    end
end

