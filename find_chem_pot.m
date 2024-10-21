function Ef = find_chem_pot(knum_tot,ek,filling,T,max_iter,threshold)

% Boltzmann constant in eV K^-1
kb = 8.617333262145e-5;

% First compute DOS at T=0
% Assumes flat bands are in the middle of the eigenspectrum
min_eval = min(ek(end/2,:));
max_eval = max(ek(end/2+1,:));

search_direction = 0;

% Number of electrons
ne = 4 + filling;

% Trial Fermi Energy
if(T==0)
   % Assume particle-hole symmetry
   Ef = 0.5*(max_eval + min_eval);
else
   Ef = 0.5*(max_eval + min_eval) + kb*T/2;
end

step = 0.5;
weight = 1/knum_tot;

% Compute Fermi level (T=0) first
for iter = 1 : max_iter
   sum_occ = 0;
   res = 0;
   for ik = 1 : knum_tot
      %for n = size(ek,1) : -1 : 1
      for n = 0 : 1
         noccs = 4*fd(ek(end/2 + n,ik),Ef,T);
         sum_occ = sum_occ + weight*noccs;
      end
   end
   integrated_ne = sum_occ;
   res = sum_occ - ne;
   rel_res = res/ne;
   if(abs(rel_res)<= threshold)
      break
   elseif(res>0)
      Ef = Ef - step;
      if(search_direction>=0)
         step = 0.5*step;
         search_direction = -1;
      end
   elseif(res<0)
      Ef = Ef + step;
      if(search_direction<=0)
         step = 0.5*step;
         search_direction = 1;
      end
   end
end
if(abs(rel_res)>=threshold)
   Ef = 0;
   disp('WARNING! In find_fermi.m maximum number of iterations has been reached, but accuracy above threshold');
   disp(join(['         for filling = ',num2str(filling)]));
end

end
