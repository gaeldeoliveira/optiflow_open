function [SD5, SD5_H] = sd5_from_h(BD5, BD5_XD5, XD5_H)

      % SUBROUTINE SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
      % IMPLICIT REAL (A-Z) 
      % ! Make dimensional shape function from bernstein polynomials
      SD5   = BD5;
      SD5_H = BD5_XD5 * XD5_H;
      
      %RETURN
      %END
      
end
      
      
% For later reference
%       SUBROUTINE SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
%       IMPLICIT REAL (A-Z) 
%       ! Make dimensional shape function from bernstein polynomials
%       SD5   = BD5
%       SD5_H = BD5_XD5 * XD5_H
%       
%       RETURN
%       END