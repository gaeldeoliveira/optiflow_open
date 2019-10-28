function delta = delta_99(hk, rt, ue_over_nue)
% Returns delta99 from hk, rt and ue_over_nue using the h1 closure relation
% (Head's shape parameter, hh1cal). Compressibility effects are ignored
% here, we are using h = hk!

    % Compute theta from rt and ue_over_nue
    theta = rt ./ ue_over_nue;

    % Compute h1 from hh1cal closure relation
    h1 = hh1cal( hk );

    % Compute delta from h1, hk and theta
    delta = theta .* (h1 - hk);

end