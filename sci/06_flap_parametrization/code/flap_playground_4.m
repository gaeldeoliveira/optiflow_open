x = linspace(0,1,100);


x_flap_hinge    = 0.8;          % Chordwise Stance Coordinate (x/c) of flap hinge
h_flap_hinge    = 0;            % Normal Stance Coordinate (z/c) of flap hinge
theta_flap      = 10;            % Flap angle in Degrees
knee_flap       = 0.02;          % Smoothing Factor around Flap (0 none, 1 a lot lot lot, 0.1 = recomended)


% Calculate Lenght and Slope of Flap
l_flap          = 1 - x_flap_hinge;                 % Lenght of Flap (over chord)
m_flap          = tan(theta_flap*pi()/180);         % Slope of Flap Chordline

% Position of Knee points
x_knee_m = x_flap_hinge - knee_flap;                % Point Before hinge
z_knee_m = 0;                                       % z (vertical(
m_knee_m = 0;                                       % Slope
x_knee_p = x_flap_hinge + knee_flap;                % Point After hinge
z_knee_p = m_flap * knee_flap;                      %

% Equivalent Trailing Edge Center Position
x_TE_eq    = 1;                                % Corresponds to sheared case value. Will be smaller in soft rotation case
z_TE_eq    = m_flap * l_flap;                  % Equivalent TE Vertical Displacement (in sheared case, real case will be smaller)

% Generate coordinate transformation definition line Interpolant Keypoints
x_pchip = [0 , x_knee_m, x_knee_p, x_TE_eq];
m_pchip = [0 , z_knee_m, z_knee_p, z_TE_eq];
z_pchip = [0 , z_knee_m, z_knee_p, z_TE_eq];
% Create a shape preserving cubic interpolant for the coordinate transformation definition line 
pp_transformation = pchip(x_pchip, z_pchip);

% Compute its first and second derivatives
d1pp_transformation = fnder(pp_transformation);
d2pp_transformation = fnder(pp_transformation, 2);

figure(1)
plot(x, ppval(pp_transformation, x));
grid on
figure(2)
plot(x, ppval(pp_transformation, x) , ...
        x, ppval(d1pp_transformation, x) , ...
            x, ppval(d2pp_transformation, x));
grid on
legend('Original', '1st Derivative' , '2nd Derivative')



%% Develop a polynomial knee function g(eta) ensuring 2nd derivative continuity between two linear branches
%   g and eta will be such that:
%       eta = (x-x1) / (x2-x1)      ,       delta = (x2-x1) / 2
%
%       f = 2 * m * delta * g(eta(x)) 
%
% And hence g will have the following boundary conditions
%
%   g(0)    = 0    ;   g(1)     = 1/2  ;
%   d1g(0)  = 0    ;   d1g(1)   = 1    ;
%   d2g(0)  = 0    ;   d2g(1)   = 0    ;
%
% The function has 6 contraints, so the polynomial must have order 6 and
% degree 5. Let us make a matrix to fit these coefficients by departing
% from the dot product form of the polynomials (in canonical basis):
%   p(x) = [1 , x , x^2 , x^3 , x^4 , x^5] * [a0 ; a1 ; a2; a3; a4; a5]
% 
% So we make six lines of these and solve for the a's!
syms x 
eta_line = [1 , x , x.^2 , x.^3 , x.^4 , x.^5];
d1eta_line = diff(eta_line, x);
d2eta_line = diff(d1eta_line, x);


% Make eta lines of matrix for the origin
eta_line_0      = subs(eta_line, 0);
d1eta_line_0    = subs(d1eta_line, 0);
d2eta_line_0    = subs(d2eta_line, 0);
% Make eta lines of matrix for eta = 1 
eta_line_1      = subs(eta_line, 1);
d1eta_line_1    = subs(d1eta_line, 1);
d2eta_line_1    = subs(d2eta_line, 1);


% Make our matrix
M = [eta_line_0 ; d1eta_line_0 ; d2eta_line_0 ; eta_line_1 ; d1eta_line_1 ; d2eta_line_1];
% We verify that it is not singular by confirming rank(M) = 6! 

% And our boundary conditions B vector. For the M matrix ordering we chose
% this yields:
%   B = [ g(0) ; d1g(0) ; d2g(0) ; g(1) ; d1g(1) ; d2g(1)];
B = [0 ; 0 ; 0 ; 1/2 ; 1 ; 0];

% Solve for the coefficients by gaussian elimination M\B or direct
% inversion inv(M)*B
a_vector = M\B;
 
% and we get a_vector = [ 0, 0, 0, 1, -1/2, 0] which is equivalent to say
% g(eta) = -1/2*x^4 + x^3  and I have verified with pen and paper that this matches all our boundary conditions! 
% So now that we checked things by hand, we can do that with the symbolic
% toolbox:
g   = eta_line * a_vector;                              % Yields : - x^4/2 + x^3
d1g = diff(  g, x);                                     % Yields : - 2*x^3 + 3*x^2
d2g = diff(d1g, x);                                     % Yields : - 6*x^2 + 6*x

% And plot this stuff
x_range = linspace(0, 1, 100);
double(subs(g, x_range))
figure(1)
plot( x_range, double(subs(g, x_range)) );
grid on
figure(2)
plot( x_range, double(subs(g, x_range)) , ...
        x_range, double(subs(d1g, x_range)) , ...
            x_range, double(subs(d2g, x_range)) );
        grid on
legend('Original', '1st Derivative' , '2nd Derivative')

%% Now calculate curvature expression analytically with the symbolic toolbox
%
% A classical formula for curvature, kappa, is given as : 
%           kappa = y'' / (1+(y')^2 )^(3/2)
% Which is not the most general form (usually parametric) that works for
% any axis orientation, but as this is a fairly well behaved case, the
% above formulat should work.
%

kappa = d2g / (1 + d1g^2)^(3/2);                % Yields: (- 6*x^2 + 6*x)/((6*x - 6*x^2)^2 + 1)^(3/2) 

% And plot curvature
figure(3)
plot( x_range, double(subs(kappa , x_range)) );
grid on


%% Now proceed to definie curvilinear coordinate transformation that preserves the orthogonality of base vectors
% Follow Modern Developments in Fluid Dynamics, Goldstein, 1938, page 119
% and define:
%   kappa - curvature
%   h1 = 1 + kappa * y      % For h1 and h2 definition see Analytical Fluid Dynamics, Emmanuel, Appendix A for a modern approach or Lamé, Leçon sur les coordonnés curvilignes, 1859, accessible on Gallica for a beautifully written classic!  asd
%   h2 = 1                  % 
%   dp = h1 * dx            % p is the transformed chordwise coordinate
%   dq  = h2 * dy           % q is the transformed normal (flapwise! ahahah!) coordinate

% Declare transformation components
syms x y 
h1 = 1 + kappa * y;
h2 = 1;

% So we can write dp
syms dx dy dp dq
dp = h1*dx;                 % Provides this: dx*((y*(- 6*x^2 + 6*x))/((6*x - 6*x^2)^2 + 1)^(3/2) + 1)
dq = h2*dy;                 % Obviously: dq = dy  <=> q = y   ... In other words the normal coordinate is invariant in this transformation!

% Things are straighforward for dq, but harder for dp. Here, we use some
% insight and do some stuff by hand:
%       dp = h1*dx      <=>     dp = 1 + k(x)*y
% Considering that for p0 = (p : x==0) and that (y orth x)
%  <=>  int(dp, p0, p) = int(1 + k(x)*y, 0, x)
%  <=>  p - p0 = x + y * int(k(x), 0, x)
% So we obtain an expression that is nice to integrate! Furthermore, in our
% particular case, because k(0) = 0, and y'(0) = 0 we get p0=0 !

% So define can now integrate kappa numerically or analytically! For this
% particular case. For this case, there exists an analytical definition of 
% the curvature integral in x, and the symbolic toolbox is able to find it 
int_kappa = int(kappa, 0, x);            % Notice that in this case int(kappa, 0, x) - int(kappa, x) = 0 !

% So we can finally define our affine transformation explicitly
p = x - int_kappa * y;                   % And here we are!!!
q = int(h2, y);                          % For this case, h2=1, we have int(h2,y)=y so we could have written just q = y; 

%% Trying our affine transformation on a grid

% Let's now generate a cartesian grid [X,Y] in the original (x,y) coordinates
x_range = linspace(0, 1, 10);
y_range = linspace(-0.1, 0.1, 10);
[X, Y] = meshgrid(x_range, y_range);

% And transform it into a grid [P, Q] in the transformed (p, q) coordinates
% Allocate memory
P = zeros(size(X));    % If the mesh is correct, then size(X) == size(Y)  
Q = zeros(size(Y));    % Same comment as line above!!!

for n_x = 1:size(X, 2)
    for n_y = 1:size(Y, 1)
        P(n_y, n_x) = double(subs(subs(p, x, X(n_y, n_x)), y, Y(n_y, n_x)));
        Q(n_y, n_x) = double(subs(g, x, X(n_y, n_x))) + double(subs(q, y, Y(n_y, n_x)));
    end
end

% And plot
figure(4)
surf(P, Q, zeros(size(P)));
hold on
plot( x_range, double(subs(g, x_range)) , 'b');
view(2)
axis([0 , 1.5, -0.75 ,0.75])
grid on
legend('Deformed Grid' , 'Deformation Line', 'Location', 'SouthEast')

% Mooij! Mooij! Mooij!


%% Let's now try to scale our transformatoion
% The knee functoin is written in the chord scale as:
%       
%  p(X) = 2*delta*m*g(eta(X)) ;     with x = eta(X) = x-x1 / x1 - x2;
%
%
%
syms m delta
