function [drhoG_vec_1,drhoG_vec_2,drhoG_vec_3,drhoG_vec_4] = compute_rho(allbands1,eigvecs,connect_m_den,tot_dim,knum_tot,CN,Ef)

%%%
%Added to split up the contributions into layers and sublattices
%%%
drhoG_vec_1(1:6) = 0;
drhoG_vec_2(1:6) = 0;
drhoG_vec_3(1:6) = 0;
drhoG_vec_4(1:6) = 0;

iG2_index = zeros(6,tot_dim);
for iq = 1 : 6
   tmp_H = kron(connect_m_den(:,:,iq),eye(2));
   tmp_H = [tmp_H,zeros(size(tmp_H));zeros(size(tmp_H)),tmp_H];
   for iG1 = 1 : tot_dim
      if(~isempty(find(tmp_H(iG1,:))))
         iG2_index(iq,iG1) = find(tmp_H(iG1,:));
      end
   end
end

for iq = 1 : 6
   for ik = 1 : knum_tot
      occ_index = find(allbands1(:,ik) < Ef);% <= Ef & allbands1(:,ik) > CN);
      
      for in = 1 : length(occ_index)
         for iG1 = 1 : tot_dim
             if(iG2_index(iq,iG1)~=0)%Is there a quicker way to do this? 
                 %also save the index where it is not zero?
                 
               if iG1 < tot_dim/2 + 1%is this correct? or do I need to add 1 or something?
                   %need ot check if iG1 is even or odd ... 
                  if rem(iG1, 2) == 0
                      %I don't quite understand why the occ_index part is
                      %the same for both G's ... 
                     drhoG_vec_2(iq) = drhoG_vec_2(iq) + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG2_index(iq,iG1),occ_index(in),ik));
                     
                     %c1 = c1 + 1;
                  else 
                     drhoG_vec_1(iq) = drhoG_vec_1(iq) + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG2_index(iq,iG1),occ_index(in),ik));
                     %c2 = c2 + 1;
                     %fprintf('iG1=%2.2f , iG2=%2.2f \n', iG1,iG2_index(iq,iG1));
                     %fprintf('iG1=%2.2f , iG2=%2.2f \n', rem(iG1, 2),rem(iG2_index(iq,iG1), 2));
                     %c1 = 0;

                  end

               else
                  if rem(iG1, 2) == 0
                     drhoG_vec_4(iq) = drhoG_vec_4(iq) + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG2_index(iq,iG1),occ_index(in),ik));
                     %c3 = c3 + 1;
                     %c2 = 1;
                  else
                     drhoG_vec_3(iq) = drhoG_vec_3(iq) + eigvecs(iG1,occ_index(in),ik)*conj(eigvecs(iG2_index(iq,iG1),occ_index(in),ik));
                     %c4 = c4 + 1;
                     %fprintf('iG1=%2.2f , iG2=%2.2f \n', iG1,iG2_index(iq,iG1));
                     %fprintf('iG1=%2.2f , iG2=%2.2f \n', rem(iG1, 2),rem(iG2_index(iq,iG1), 2));
                     %c2 = 0;
                  end
                  
               end 
            end
         end
      end
   end
end
end
