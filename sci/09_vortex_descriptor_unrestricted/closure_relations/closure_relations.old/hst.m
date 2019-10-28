function [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq)
%HST Turbulen HS correlation
%       new correlation  29 Nov 91
%       (from  arctan(y+) + Schlichting  profiles)
% Note: HSTIVW is currently not used by Rfoilsuc
[ hs, hs_hk, hs_rt, hs_msq  ] = hstivw2( hk, rt, msq);

end