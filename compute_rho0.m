function [drho0_1,drho0_2,drho0_3,drho0_4] = compute_rho0(allbands1,eigvecs,tot_dim,knum_tot,CN,Ef)

drho0_1 = 0;
drho0_2 = 0;
drho0_3 = 0;
drho0_4 = 0;

for ik = 1 : knum_tot
  % Occupied bands only     
  occ_index = find(allbands1(:,ik) < Ef);% & allbands1(:,ik) > CN);

  for in = 1 : length(occ_index)
     for iG1 = 1 : tot_dim
           if iG1 < tot_dim/2+1%is this correct? or do I need to add 1 or something?
              if rem(iG1, 2) == 0
                 drho0_2 = drho0_2 + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG1,occ_index(in),ik));
              else
                 drho0_1 = drho0_1 + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG1,occ_index(in),ik));
              end
           else
              if rem(iG1, 2) == 0
                 drho0_4 = drho0_4 + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG1,occ_index(in),ik));
              else
                 drho0_3 = drho0_3 + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG1,occ_index(in),ik));
              end
           end 
     end
  end
end

end
