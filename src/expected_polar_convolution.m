function [alpha_bar_range, ld_expected] = ...
    expected_polar_convolution(sigma, n_sigma, n_alpha_til, alpha_range, ld_alpha)
% Typical inputs
% sigma       = 4;              % [deg  ] - Standard deviation of angle of attack perturbations
% n_alpha_til =  64;            % [int. ] - number of elements convolution integral kernel
% n_sigma     = 1.64485;        % [adim.] - Width of integration range (in sigmas, affects pdf/cdf coverage)
% alpha_range = ap.alpha_range  % [deg  ] - range over which aero data is defined from steady state simulation
% ld_alpha    = experiment.ap.cl_alpha(alpha_range) ./ experiment.ap.cd_alpha(alpha_range);

% % Check that sigma is not zero
if sigma > 0
    % If it is not, then let's try to convolute
    
    % % Inform about pdf converage
    % pdf_coverage = normcdf(n_sigma*sigma ,0,sigma) - normcdf(-n_sigma*sigma,0,sigma); disp(['PDF coverage = '  num2str(pdf_coverage)]);
    
    % % Make a function for L over D at demanded angle of attack
    ld_fun   = @(alpha) interp1(alpha_range, ld_alpha, alpha);
    % Determine range of mean angles of attack around which convolution is
    % doable (function evaluable over entire range)
    alpha_bar_min = min(alpha_range) + sigma * n_sigma;
    alpha_bar_max = max(alpha_range) - sigma * n_sigma;
    % Make range of mean angles of attack for convolution
    % alpha_bar_range = linspace(alpha_bar_min, alpha_bar_max, n_alpha_bar);
    alpha_bar_range = alpha_range(and(alpha_range >= alpha_bar_min, alpha_range <= alpha_bar_max));
    % Allocate space for storing results of convolution
    ld_expected     = zeros(size(alpha_bar_range));
    
    % % Check that everything is doable
    if alpha_bar_max > alpha_bar_min
        % If doable, go on!
        % % Generate probability kernel
        %   alpha_til  - instantaneous angle of attack
        %   sigma      - standard deviation of instantaneous angle of attack perturbations
        p_alpha_til_fun    = @(alpha_til, sigma) normpdf(alpha_til, 0, sigma);
        % Generate centered range for convolution kernel
        alpha_til_kernel   = linspace(-n_sigma*sigma, n_sigma*sigma, n_alpha_til);
        % Generate centered pdf   for convolution kernel
        p_alpha_til_kernel = p_alpha_til_fun(alpha_til_kernel, sigma);
        % Estimate kernel loss
        p_alpha_til_kernel_loss = trapz(alpha_til_kernel, p_alpha_til_kernel);
        % Normalize kernel to limit loss
        p_alpha_til_kernel_normalized = p_alpha_til_kernel /p_alpha_til_kernel_loss;
        
        % % Move on to compute convolution
        for n_alpha_bar = 1:length(alpha_bar_range)
            % Get current mean angle of attack
            alpha_bar                = alpha_bar_range(n_alpha_bar);
            % Make shifted kernel range for convolution
            alpha_til_conv           = alpha_bar + alpha_til_kernel;
            % Make l over d at kernel positions
            f_conv                   = ld_fun(alpha_til_conv);
            % Compute product of l over d with probability kernel
            f_p_conv                 = f_conv .* p_alpha_til_kernel_normalized;
            % Integrate to obtain expected value of l over d
            ld_expected(n_alpha_bar) = trapz(alpha_til_kernel, f_p_conv);
        end
        
    else
        % Not doable, notify it!
        alpha_bar_range = [0 1];
        ld_expected = [NaN NaN];
    end
    
else
    % Sigma was 0, expected polar corresponds to original one!
    alpha_bar_range = alpha_range;
    ld_expected     = ld_alpha;
end

end
