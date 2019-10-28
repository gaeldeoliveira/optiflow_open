function [ di, di_hs, di_us, di_cf, di_st ] = dit( hs, us, cf, ctau)
%DIT    Turbulent dissipation function  ( 2 CD/H* )
%   From Rfoilsuc

di    =   ( 0.5 * cf .* us + ctau .* (1.0-us) ) * 2.0 ./hs;
di_hs = - ( 0.5 * cf .* us + ctau .* (1.0-us) ) * 2.0 ./hs^2;
di_us =   ( 0.5 * cf       - ctau             ) * 2.0 ./hs;
di_cf =   ( 0.5   *     us                     ) * 2.0 ./hs;
di_st =   (                       .* (1.0-us) ) * 2.0 ./hs;


end


% Original Rfoilsuc code
%       SUBROUTINE DIT( HS, US, CF, ST, DI, DI_HS, DI_US, DI_CF, DI_ST )
% !
% !---- Turbulent dissipation function  ( 2 CD/H* )
%       DI    =  ( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS
%       DI_HS = -( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS**2
%       DI_US =  ( 0.5*CF    - ST*ST          ) * 2.0/HS
%       DI_CF =  ( 0.5   *US                  ) * 2.0/HS
%       DI_ST =  (            2.0*ST*(1.0-US) ) * 2.0/HS
% !
%       RETURN
%       END
