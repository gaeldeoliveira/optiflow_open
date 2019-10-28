function [BD5, BD5_X] = bst_d5(X, A1, A2, A3, A4, A5, A6)
      % SUBROUTINE BST_D5( X, A1, A2, A3, A4, A5, A6, BD5, BD5_X)
      % IMPLICIT REAL (A-Z) 
      
      % ! Make Bernstein polynomial basis for degree 5 (order 6)
      BD5_R0 =         X.^5              ;
      BD5_R1 =  5.0 .* X.^4 .* (1.0-X)   ;
      BD5_R2 = 10.0 .* X.^3 .* (1.0-X).^2;
      BD5_R3 = 10.0 .* X.^2 .* (1.0-X).^3;
      BD5_R4 =  5.0 .* X    .* (1.0-X).^4;
      BD5_R5 =                 (1.0-X).^5;
      % ! Combine basis linearly
      BD5    = A1.*BD5_R0 + A2.*BD5_R1 + A3.*BD5_R2 + A4.*BD5_R3 + A5.*BD5_R4 + A6.*BD5_R5;
      
      % ! Now make derivatives of polynomials from Bernstein basis of degree 5 (order 6)
      BD5_R0_X =  1.0*(5.0 *X.^4                                        );
      BD5_R1_X =  5.0*(4.0 *X.^3 .*(1.0-X)    -  1.0 .*X.^4             );
      BD5_R2_X = 10.0*(3.0 *X.^2 .*(1.0-X).^2 -  2.0 .*X.^3 .*(1.0-X)   );
      BD5_R3_X = 10.0*(2.0 *X    .*(1.0-X).^3 -  3.0 .*X.^2 .*(1.0-X).^2);
      BD5_R4_X =  5.0*(1.0       .*(1.0-X).^4 -  4.0 .*X    .*(1.0-X).^3);
      BD5_R5_X =  1.0*(                       -  5.0        .*(1.0-X).^4);
      % ! Combine derivatives of basis linearly
      BD5_X    = A1.*BD5_R0_X + A2.*BD5_R1_X + A3.*BD5_R2_X + A4.*BD5_R3_X + A5.*BD5_R4_X + A6.*BD5_R5_X;
      
%!      BD5_R0_X = 
      
%!      BD5   = A1*X +     A2*X**2
%!      BD5_X = A1   + 2.0*A2*X
      
%     RETURN
%     END

end



% For reference
% SUBROUTINE BST_D5( X, A1, A2, A3, A4, A5, A6, BD5, BD5_X)
%       IMPLICIT REAL (A-Z) 
%       
%       ! Make Bernstein polynomial basis for degree 5 (order 6)
%       BD5_R0 =        X**5
%       BD5_R1 =  5.0 * X**4 * (1.0-X)
%       BD5_R2 = 10.0 * X**3 * (1.0-X)**2
%       BD5_R3 = 10.0 * X**2 * (1.0-X)**3
%       BD5_R4 =  5.0 * X    * (1.0-X)**4
%       BD5_R5 =               (1.0-X)**5
%       ! Combine basis linearly
%       BD5    = A1*BD5_R0 + A2*BD5_R1 + A3*BD5_R2 + A4*BD5_R3 + A5*BD5_R4 + A6*BD5_R5
%       
%       ! Now make derivatives of polynomials from Bernstein basis of degree 5 (order 6)
%       BD5_R0_X =  1.0*(5.0 *X**4                                     )
%       BD5_R1_X =  5.0*(4.0 *X**3 *(1.0-X)    -  1.0 *X**4            )
%       BD5_R2_X = 10.0*(3.0 *X**2 *(1.0-X)**2 -  2.0 *X**3 *(1.0-X)   )
%       BD5_R3_X = 10.0*(2.0 *X    *(1.0-X)**3 -  3.0 *X**2 *(1.0-X)**2)
%       BD5_R4_X =  5.0*(1.0       *(1.0-X)**4 -  4.0 *X    *(1.0-X)**3)
%       BD5_R5_X =  1.0*(                      -  5.0       *(1.0-X)**4)
%       ! Combine derivatives of basis linearly
%       BD5_X    = A1*BD5_R0_X + A2*BD5_R1_X + A3*BD5_R2_X + A4*BD5_R3_X + A5*BD5_R4_X + A6*BD5_R5_X
%       
% !      BD5_R0_X = 
%       
% !      BD5   = A1*X +     A2*X**2
% !      BD5_X = A1   + 2.0*A2*X
%       
%       RETURN
%       END