function [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq)
% Wrapper for pseudo vectorized execution of hst_seq closure (HST
% subroutine in Rfoil Fortran Code)
    hs     = zeros(size(hk));
    hs_hk  = zeros(size(hk));
    hs_rt  = zeros(size(hk));
    hs_msq = zeros(size(hk));
    for n_hk = 1:length(hk)
       [ hs(n_hk), hs_hk(n_hk), hs_rt(n_hk), hs_msq(n_hk)] = hst_seq( hk(n_hk), rt, msq);
    end
end