function [ di, di_hs, di_us, di_cf, di_st , di_sv] = dits( hs, us, cf, st, sv)
%DITS 
%       Gael de Oliveira : Modify Dissipation Coefficient to account for vertical velocity term (suction) in shear
%                    This is the modification of Merchant

% !---- Turbulent dissipation function  ( 2 CD/H* )
      di    =  ( 0.5 * cf .* us  + 0.5 * sv .* us.^2 + (st.^2) .* (1 - us) ) .* 2.0 ./ hs;
% Derivatives      
      di_hs = -( 0.5 * cf .* us  + 0.5 * sv .* us.^2 + (st.^2) .* (1 - us) ) .* 2.0 ./ (hs.^2);
      di_us =  ( 0.5 * cf        + 1.0 * sv .* us    - (st.^2)             ) .* 2.0 ./ hs;
      di_cf =  ( 0.5       * us                                            ) .* 2.0 ./ hs;
      di_st =  (                                        2 * st .* (1 - us) ) .* 2.0 ./ hs;
      di_sv =  (                 + 0.5      .* US.^2                       ) .* 2.0 ./ hs;

end


% %% Original Fortran Code
% ! Gael de Oliveira : Modify Dissipation Coefficient to account for vertical velocity term (suction) in shear
% !                    This is the modification of Merchant
% ! 	  Original DIT subroutine
% 
% ! 	  New DITS subroutine, accounting for suction :
%       SUBROUTINE DITS( HS, US, CF, ST, SV, DI, DI_HS, DI_US, DI_CF, DI_ST , DI_SV)
% !
% !---- Turbulent dissipation function  ( 2 CD/H* )
%       DI    =  ( 0.5*CF*US  + 0.5 * SV * US**2 + ST*ST*(1.0-US) ) * 2.0/HS
%       DI0    =  ( 0.5*CF*US + 0.5 * 0 * US**2 + ST*ST*(1.0-US) ) * 2.0/HS
% !      WRITE(*,*) 'DI = ' , DI , '  DI0 =' , DI0 , '   SV =' , SV
%       DI_HS = -( 0.5*CF*US  + 0.5 * SV * US**2 + ST*ST*(1.0-US) ) * 2.0/HS**2
%       DI_US =  ( 0.5*CF     + 1.0 * SV * US    - ST*ST          ) * 2.0/HS
%       DI_CF =  ( 0.5   *US                                     ) * 2.0/HS
%       DI_ST =  (                               2.0*ST*(1.0-US) ) * 2.0/HS
%       DI_SV =  (            + 0.5      * US**2                  ) * 2.0/HS
% !
%       RETURN
%       END
% 
