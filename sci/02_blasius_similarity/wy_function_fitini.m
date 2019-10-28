        function wy = wy_function_fitini(y, b)
            % wY weighting function, as described in report.
            % Supports vectorized evaluation (conditionality handled as
            % multiplication!)
            % Taken from the plasma_descriptor class to be nice with
            % Fitini, la gentille équation!
            wy = 0.5*pi() * sin( 0.5*pi()*(y/b + 1)) .* (y/b<1) .* (y>0);
        end