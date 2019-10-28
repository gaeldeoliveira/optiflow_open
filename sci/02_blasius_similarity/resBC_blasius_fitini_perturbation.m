
function res = resBC_blasius_fitini_perturbation(ya,yb)
% Boundary conditions.
% As described in the report as Rbc
res = [  ya(1)
    ya(2)
    yb(2)];
end