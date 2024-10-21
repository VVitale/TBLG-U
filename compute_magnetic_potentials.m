function [V0p,V1p] = compute_magnetic_potentials(U_param,potname,knum_tot,up0_1,up0_2,up0_3,up0_4,dw0_1,dw0_2,dw0_3,dw0_4,...
    up_1,up_2,up_3,up_4,dw_1,dw_2,dw_3,dw_4);

if(strcmp(potname,'NONE'))
    V0p = 0;
    V1p = 0;
elseif(strcmp(potname,'MAFM'))
    V_UP0p = U_param*2*real(sum(up0_1)+sum(up0_3)+sum(dw0_2)+sum(dw0_4))/6/4/knum_tot;
    V_DW0p = U_param*2*real(sum(dw0_1)+sum(dw0_3)+sum(up0_2)+sum(up0_4))/6/4/knum_tot;

    V_UPp = U_param*2*real(sum(up_1)+sum(up_3)+sum(dw_2)+sum(dw_4))/6/4/knum_tot;
    V_DWp = U_param*2*real(sum(dw_1)+sum(dw_3)+sum(up_2)+sum(up_4))/6/4/knum_tot;

    V0p = abs(V_UP0p - V_DW0p);
    V1p = abs(V_UPp - V_DWp);

elseif(strcmp(potname,'NAFM'))
    V_UPp = U_param*2*real(sum(up_1)+sum(up_3)+sum(dw_2)+sum(dw_4))/6/4/knum_tot;
    V_DWp = U_param*2*real(sum(dw_1)+sum(dw_3)+sum(up_2)+sum(up_4))/6/4/knum_tot;

    V0p = 0;
    V1p = abs(V_UPp - V_DWp)

elseif(strcmp(potname,'HAFM'))%Need to check this one ... 
    V_UPp = U_param*2*imag(sum(up_1)+sum(up_3)+sum(dw_2)+sum(dw_4))/6/4/knum_tot;
    V_DWp = U_param*2*imag(sum(dw_1)+sum(dw_3)+sum(up_2)+sum(up_4))/6/4/knum_tot;

    V0p = 0;
    V1p = abs(V_UPp - V_DWp);

elseif(strcmp(potname,'FM'))
    V_UPp = U_param*2*real(sum(up_1)+sum(up_2)+sum(up_2)+sum(up_4))/6/4/knum_tot;
    V_DWp = U_param*2*real(sum(dw_1)+sum(dw_3)+sum(dw_2)+sum(dw_4))/6/4/knum_tot;

    V0p = 0;
    V1p = abs(V_UPp - V_DWp);

end
end
