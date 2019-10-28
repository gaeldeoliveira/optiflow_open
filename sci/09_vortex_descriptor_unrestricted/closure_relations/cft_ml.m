function [cf, cf_hk, cf_rt , cf_msq] = cft_ml( hk, rt, msq, SPR)
% Pseudo-vectorization wrapper for machine learning skin friction relation
% (clone of fortran Rfoil implementation

        cf     = zeros(size(hk));
        cf_hk  = zeros(size(hk));
        cf_rt  = zeros(size(hk));
        cf_msq  = zeros(size(hk));
        for n_hk = 1:length(hk)
            [cf(n_hk), cf_hk(n_hk), cf_rt(n_hk), cf_msq(n_hk)] = cft_ml_seq( hk(n_hk), rt, msq, SPR);
        end

end

