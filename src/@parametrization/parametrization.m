classdef parametrization < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        order
        C_t
        degree
    end
    
    properties (SetAccess = private) 
      % Vectors
      s_vector
      sc_vector
      sc_vector_diff 

      % Matrices for fitting through use of the 'order-2' inner points 
      tx_positions
      s_matrix
      sc_matrix
      
      % Non singular (invertible) part of SC matrix and RHS terms
      sc_submatrix
      sc_submatrix_inverse
      sc_a0_line
      sc_an_line
      tx_positions_inner
      fitstring = ''
      
      M1               % horizontal vector of Bernstein Coefficients * M = vertical vector of standard polynomial coefficients (for the shape function)
      % M_inverse      % Minv * vertical vector of standard polynomial coefficients = vertical vector of Bernstein Coefficients (for the shape function)
      D                % horizontal vector of standard polynomial coefficients * D = horizontal vector of polynomial coefficients of derivative polynomial
      
    end
    
    methods      
        function p = parametrization(name , order, varargin)
            p.order = order;
            p.degree = order - 1;
            p.name = name;            
            
            if size(varargin,2) > 0
               syms t;
               p.C_t = varargin{1};
            else
               syms t;
               p.C_t = (1-t).*sqrt(t);
            end
            
            p.tx_positions = 0:(1/p.degree):1;
            p.tx_positions_inner = p.tx_positions(2:p.degree);
            
            p.make_vectors();
            p.make_matrices();
            p.extract_non_singular_system();
            p.invert_non_singular_matrix();
            
            p.make_fitstring_primitive();
            
            p.make_standard_form_conversion_matrix();
            p.make_polynomial_derivation_matrix();
            
        end        
        
        function b = Binomial(p,i)
            N = p.degree;
            b = factorial(N)/(factorial(N-i)*factorial(i));
        end      
        function s = S1(p , t, i)
            N1 = p.degree;
            s = p.Binomial(i) *(t.^i).*(1-t).^(N1-i);
        end        
        function make_vectors(p)
            syms t;
            s_vector = [];
            
            for k = 0:p.degree
                s_vector = [ s_vector , p.S1(t, k) ] ;
            end
            
            p.s_vector = s_vector;
            p.sc_vector = p.C_t * s_vector;
            p.sc_vector_diff =  diff(p.sc_vector,t);
        end        
        function s_line = s_line(p, t)
            s_line = subs(p.s_vector, 't' , t);
        end
        function sc_line = sc_line(p, t)
            sc_line = subs(p.sc_vector, 't' , t);
        end  
        function make_matrices(p)            
            
            s_matrix = [];
            sc_matrix = [];
            
            syms t;
            
            for k = (0:p.degree)+1
                s_matrix  = [s_matrix  ; p.s_line(p.tx_positions(k))];
                sc_matrix = [sc_matrix ; p.sc_line(p.tx_positions(k))];
            end
            
            p.s_matrix =  s_matrix;
            p.sc_matrix = sc_matrix;
            
        end
        function extract_non_singular_system(p)
            N1 = p.degree;
            
            p.sc_submatrix = p.sc_matrix(2:N1, 2:N1);
            p.sc_a0_line =   p.sc_matrix(2:N1, 1);
            p.sc_an_line =   p.sc_matrix(2:N1, N1+1);
            
        end
        function invert_non_singular_matrix(p)
            p.sc_submatrix_inverse = p.sc_submatrix^(-1);
        end
        function make_fitstring_primitive(p)
            % Make the string primitive for the fitting function fo
            % non-linear least square fitting of middle coefficients
            
            p.fitstring = [ '*(' ,  char(p.sc_vector(1)) , ')'];
            
            for n_poly = 2:(length(p.sc_vector)-1)
                parcel = ['a', num2str(n_poly) ,'*(' ,  char(p.sc_vector(n_poly)) , ')'];
                p.fitstring = [p.fitstring , ' + ' , parcel ];
            end
            
            parcel = [ ' (' , char(p.sc_vector(end)) , ')*'];
            p.fitstring = [p.fitstring , ' + ' , parcel];
        end
        function ffstring = make_fitfunction(p, a0, aend)
            ffstring = [num2str(a0) , p.fitstring , num2str(aend)];
        end
        function ffstring = make_fitfunction_full_non_linear(p)
            ffstring = ['a1' , p.fitstring , 'a', num2str(p.order)];
        end        
        function f = make_fit_object_full_non_linear(p, cst_parameters_vector)
            
            ffstring = p.make_fitfunction_full_non_linear();
            
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',-3,...
                'Upper',3, ...
                'TolX' , 1e-8, ...
                'TolFun' , 1e-8, ...
                'DiffMinChange' , 1e-10, ...
                'Startpoint',cst_parameters_vector);
            
            f = fittype(ffstring , 'independent', {'t'} , 'options' , s);
            % Comme un Manouche Sans Guitare - Thomas Dutronc
        end
        function f = make_fit_object(p, cst_parameters_vector)
            ffstring = p.make_fitfunction(cst_parameters_vector(1), cst_parameters_vector(end));
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',-3,...
                'Upper',3, ...
                'TolX' , 1e-8, ...
                'TolFun' , 1e-8, ...
                'DiffMinChange' , 1e-10, ...
                'Startpoint',cst_parameters_vector(2:(end-1)));  % startpoint is one on every
            f = fittype(ffstring , 'independent', {'t'} , 'options' , s);
        end
        function make_standard_form_conversion_matrix(p)
            % Create symbolic independent variable
            syms t
            
            % Create transposed Matrix in which to store coefficients
            p.M1 = zeros(p.order);
            c = cell(p.order, 1);
                        
            for n = 1:p.order
                
                % Extract polynomial coefficients from n^th element of
                % s_vector (corresponding to the (n-1)^th elementary shape
                % function 
                c{n} = sym2poly(p.s_vector(n));
                
                
                % Below, an obsolete version with the wrong ordering
%                     % Extract polynomial coefficients from n^th element of s_vector
%                     c{n} = coeffs(p.s_vector(n), t);
%                     % Now padd the c{n} polynomial coefficient vector for the lower order terms
%                     c{n} = [c{n} , zeros(1 , n-1)];
%                     % Now write c{n} as a line in Mt matrix
                p.M1(n , 1:p.order) = c{n};
            end
            
%             % We actually just composed the transposed matrix (it was more convenient),
%             % so we now get the matrix we are looking for by tranposing Mt
%             
%             p.M = transpose(Mt);
%             
%             % Generate inverse for inverse operation
%             p.M_inverse = (p.M)^-1;
                        
        end
        
        function make_polynomial_derivation_matrix(p)
            p.D = zeros(p.order , p.order);
            for i = 2:p.order
                p.D(i,i-1) = p.order - (i-1);
            end
        end
        
%         function make_CST_to_PARSEC_conversion_matrix(p)
%             % Solve dimension problem here
%             
%             % Create transposed Matrix in which to store coefficients
%             Mt = zeros(p.order);
%                         
%             for n = 1:p.order
%                 % Extract polynomial coefficients from n^th element of s_vector
%                 c{n} = coeffs(p.s_vector(1)*(1-t), t);
%                 % Now padd the c{n} polynomial coefficient vector for the lower order terms
%                 c{n} = [c{n} zeros(n-1)];
%                 % Now write c{n} as a line in Mt matrix
%                 Mt(n , 1:p.order) = c{n};
%             end
%             
%             % We actually just composed the transposed matrix (it was more convenient),
%             % so we now get the matrix we are looking for by tranposing Mt
%             
%             p.M = transpose(Mt);
%             
%             % Generate inverse for inverse operation
%             p.M_inverse = (p.M)^-1;
%                         
%         end
    end
end