function [ hs, hs_hk, hs_rt, hs_msq  ] = hstivw( hk, rt, msq)
%---- Turbulent HS correlation (ENGINEERING method from R.v.Rooy)
%      fit from Thomas, comparable with Drela, DHSINF= 0.08, no Re_theta

A0 =  2.0;
A1 = -2.117532;
A2 =  2.63649;
A3 = -.824744;
A4 =  .130206;
A5 =  .015373;
A6 =  .074399;
AA =  .4342945;
HKLN = AA*log(hk);
%
hs = (A0 + A1*HKLN    + A2*HKLN.^2 + A3*HKLN.^3 ...
         + A4*HKLN.^4 + A5*HKLN.^5 + A6*HKLN.^6);
     
hs_hk = (A1*AA./hk            + A2*AA*2*HKLN./hk    +A3*AA*3*HKLN.^2./hk ...
        +A4*AA*4*HKLN.^3./hk + A5*AA*5*HKLN.^4./hk ...
        +A6*AA*6*HKLN.^5./hk);

hs_rt = hs * 0;
%---- Whitfield's minor additional compressibility correction
      FM     = 1.0 + 0.014 *msq;
      hs     = ( hs + 0.028*msq )    ./ FM;
      hs_hk  = ( hs_hk          )    ./ FM;
      hs_rt  = ( hs_rt          )    ./ FM;
      hs_msq =  0.028./FM - 0.014*hs ./ FM;


end


% %%
%        SUBROUTINE HSTIVW ( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
%       IMPLICIT REAL (A-H,M,O-Z)
% !
% !---- Turbulent HS correlation (ENGINEERING method from R.v.Rooy)
% !      fit from Thomas, comparable with Drela, DHSINF= 0.08, no Re_theta
% !
% !      Dummy variables
%       HO    = 3.0
% ! datacode -->
% ! Supress dummy argument warning in dirty way to make stuff readable
% !      HO_RT = 0.0
%       HO_RT = 0.0 * RT
% ! <-- datacode
% !
% !----- attached branch and separated branch
% !
%       A0 =  2.0
%       A1 = -2.117532
%       A2 =  2.63649
%       A3 = -.824744
%       A4 =  .130206
%       A5 =  .015373
%       A6 =  .074399
%       AA =  .4342945
%       HKLN = AA*ALOG(HK)
% !
%       HS = (A0+A1*HKLN+A2*HKLN**2+A3*HKLN**3 &
%              +A4*HKLN**4+A5*HKLN**5+A6*HKLN**6)
%       HS_HK = (A1*AA/HK+A2*AA*2*HKLN/HK+A3*AA*3*HKLN**2/HK &
%              +A4*AA*4*HKLN**3/HK+A5*AA*5*HKLN**4/HK &
%              +A6*AA*6*HKLN**5/HK)
% 
% !---- Whitfield's minor additional compressibility correction
%       FM = 1.0 + 0.014*MSQ
%       HS     = ( HS + 0.028*MSQ ) / FM
%       HS_HK  = ( HS_HK          ) / FM
%       HS_RT  = ( HS_RT          ) / FM
%       HS_MSQ = 0.028/FM  -  0.014*HS/FM
% !
%       RETURN
%       END