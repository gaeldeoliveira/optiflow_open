function [ us2, us2_h2, us2_hk2, us2_hs2] = usg(h2, hk2, hs2)
%USG Returns slip speed (us) from G-beta relation
%   See Gael's master thesis for explanations

us2     = 0.5 * hs2 .* ( 3.0 - 4.0*(hk2-1.0)./h2    )/3.0;
us2_hs2 = 0.5 *        ( 3.0 - 4.0*(hk2-1.0)./h2    )/3.0;
us2_hk2 = 0.5 * hs2 .* (     - 4.0          ./h2    )/3.0;
us2_h2  = 0.5 * hs2 .* (       4.0*(hk2-1.0)./h2.^2 )/3.0;

end

% Original Fortran Code
%       SUBROUTINE USG(H2 , HK2 , HS2 , US2, US2_H2 , US2_HK2 , US2_HS2)
%       IMPLICIT REAL (A-Z)
% !     ---- Gael de Oliveira 	(Rededuction of Drela for generic B)
% 
%       US2     = 0.5*HS2*( 3.0 - 4.0*(HK2-1.0)/H2   )/3.0
%       US2_HS2 = 0.5  *  ( 3.0 - 4.0*(HK2-1.0)/H2   )/3.0
%       US2_HK2 = 0.5*HS2*(     - 4.0          /H2   )/3.0
%       US2_H2  = 0.5*HS2*(       4.0*(HK2-1.0)/H2**2)/3.0
%       RETURN
%       END
