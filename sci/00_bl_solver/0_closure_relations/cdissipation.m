function [ cdi, cdi_us, cdi_cf, cdi_ctau ] = cdissipation( us, cf, ctau)
%DIT    Turbulent dissipation Coefficient
%   From Rfoilsuc, with modifications

cdi    =   ( 0.5 * cf .* us + ctau .* (1.0-us) );
cdi_us =   ( 0.5 * cf       - ctau             );
cdi_cf =   ( 0.5   *     us                    );
cdi_ctau =   (                        (1.0-us) );


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
