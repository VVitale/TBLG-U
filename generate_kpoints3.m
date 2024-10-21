function [all_kpts1,scale_axis1,knum_tot,vkp,qvecs] = generate_kpoints3(gstar_b1,...
    gstar_b2,chi,knum,full_bz,bz_n,plot_DOS)
    % K G M K
    kk1a=chi*(gstar_b1 + 2*gstar_b2)/3;
    kk1b=chi*(2*gstar_b1 + gstar_b2)/3;
    kk1c=chi*(-gstar_b1 + gstar_b2)/3;
    
    % Ks'
    kk2a=-chi*(gstar_b1-gstar_b2)/3;
    kk2b=-chi*(gstar_b1+2*gstar_b2)/3;
    kk2c=-chi*(-2*gstar_b1-gstar_b2)/3;
    
    kk3=0*gstar_b1;%\Gamma point
    kk4=chi*(gstar_b1+gstar_b2)/2;%and probably the M point
    %scan_klist1= chi*[kk3;kk1a;kk4;kk3];
    scan_klist1= chi*[kk1a;kk3;kk4;kk1b];
    
    scan_klist1(:,3)=0;
    
    all_kpts1 = zeros(3,knum);
    scale_axis1 = [];
    [ all_kpts1, scale_axis1] = generate_k_line( knum, scan_klist1 );
    kp = 0;
    vkp = 0;
    ikp = 0;
    qvecs = 0;
    if (full_bz == 1 || plot_DOS )
       recL = [gstar_b1;gstar_b2;[0,0]];%
       K = zeros(bz_n*bz_n,3);
       all_kpts1 = zeros(bz_n*bz_n,3);
       kxi=0;
       kyi=0;
       kxf=(bz_n-1)/bz_n;
       kyf=kxf;
       kx=linspace(kxi,kxf,bz_n);
       ky=linspace(kyi,kyf,bz_n);
       ikk = 0;
       %shift = rand(1)*norm(gstar_b1)/(10*bz_n);
       shift = 0.0;
       for ikx = 1 : bz_n
          for iky = 1 : bz_n
              ikk = ikk + 1;
              K(ikk,:) = [kx(ikx)+shift,ky(iky)+shift,0.0];
          end
       end
       scale_axis = sqrt(sum(K.^2,2));
       all_kpts1 = K*recL;
       knum_tot = bz_n*bz_n;
       clear K;
    end
    knum_tot = size(all_kpts1,1);
end
