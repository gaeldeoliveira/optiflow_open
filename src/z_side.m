function z = z_side(tx, A, z_te , C_t)
% Return the height of the parametrized profile usingthe CST method. 
% Receives a set of weighting parameters A (a vector, and not a matrix as
% in the wing)
% tx is the local coordinate varying from 0 to 1 within the profile
% (leading edge is at 0)


Ncoef1 = length(A);  % Get size of matrix 

N1 = Ncoef1 - 1; % Degree of polynomial is equal to number of coefficients minus one

Binomial = @(N, i) factorial(N)/(factorial(N-i)*factorial(i));
S1 = @(t , i) Binomial(N1, i) *(t.^i).*(1-t).^(N1-i);

C= @(tx) subs(C_t, 't',tx);  

Sn =0;

for i = 0:N1
        Sn = Sn + A(i+1).*S1(tx,i);
end

z =   double(C(tx) .* Sn + z_te * tx);
%zextra= - C(tx) .* Snextra;