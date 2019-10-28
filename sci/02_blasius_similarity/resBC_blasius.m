
function res = resBC_blasius(ya,yb)
% Boundary conditions.
% As described in the report as Rbc
res = [  ya(1)
    ya(2)
    yb(2)-1];
end