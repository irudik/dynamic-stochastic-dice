% Smolyak_Polynomial.m is a routine that constructs the multidimensional 
% basis functions of Smolyak polynomial of the approximation level that  
% corresponds to the previously constructed Smolyak (multidimensional) grid 
% points; see "Smolyak method for solving dynamic economic models: Lagrange 
% interpolation, anisotropic grid and adaptive domain" by Kenneth L. Judd, 
% Lilia Maliar, Serguei Maliar and Rafael, (2014). Journal of Economic 
% Dynamics and Control 44, 92–123 (henceforth, JMMV (2014)). 
%
% This version: August 11, 2016 edited by Ivan Rudik. First version: May 30, 2011.
% Rudik made three changes:
% 1) Pre-allocating Smol_bases
% 2) Store vectors in pre-allocated Smol_bases matrix instead of
% concatenation
% 3) Allow for computing derivatives
% These two changes allow for the evaluation of future state in value
% function iteration to be sped up several-fold.
% -------------------------------------------------------------------------
% Inputs:  "points"    is the matrix of points in which the polynomial basis 
%                      functions must be evaluated; numb_pts-by-d
%          "d"         is the number of dimensions (state variables)
%          "Smol_elem" is the matrix of the subindices of the Smolyak
%                      unidimensional elements; these elements can be either 
%                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
%                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
%                      in the former case, the indices i1,...,id that jointly 
%                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
%                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
%                      in the later case, they are a subset of the above 
%                      indices 
%
% Output:  "Smol_bases" is the matrix of multidimensional basis functions of 
%                      Smolyak polynomial of the given level of approximation, 
%                      evaluated in data matrix "points"
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


 function Smol_bases = Smolyak_Polynomial(points,d,mu,Smol_elem,deriv,dist)
 if nargin < 5
     deriv = zeros(1,6);
 end
% Smolyak polynomial is given by the sum of multidimensional basis functions, 
% multiplied by the coefficients; see formula (15) in JMMV (2014). By 
% convention, the first basis function is given by 1 (unity). 

% Unidimensional basis functions are given by Chebyshev polynomial bases; 
% in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
% degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
% this notation here

% 1. Construct the unidimensional basis functions and evaluate them in
% all the points of matrix "points"
% -------------------------------------------------------------------------
i_max = mu+1;                   % The maximum subindex of unidimensional
                                % set S_i whose points are used to
                                % construct Smolyak grid; e.g., for mu=1, 
                                % we consider up to S_i_max={0,-1,1} where
                                % i_max=mu+1=1+1=2
                                
% Compute the number of elements in the i_max-th unidimensional set of 
% elements S_i_max; this coincides with a maximum subindex of elements 
% (unidimensional grid point or unidimensional basis function); 
% see Section 2.2.1 in JMMV (2014)

m_i_max(i_max==1) = 1;          % If i_max=1, then m(i_max)=1, i.e., set  
                                % S_1={0} and the maximum subindex is 1
m_i_max(i_max>1) =  2.^(i_max(i_max>1)-1) + 1;
                                % If i_max>1, then m(i_max)= 2^(i_max-1)+1;
                                % e.g., for S_2={0,-1,1}, the maximum
                                % subindex is 3

                                 
numb_pts = size(points,1);      % Compute the number of points (rows),   
                                % "numb_pts" in the matrix of points 
                                % "points", in which the polynomial bases  
                                % must be evaluated              
phi = ones(numb_pts,d,m_i_max); 
                                % Allocate memory to a matrix of the
                                % unidimensional polynomial bases "phi_n",  
                                % evaluated in all the points 
                     
% For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
% our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
% construction
                               
phi(:,:,2) = points;            % For a polynomial bases "phi_n" with n=2, 
                                % we have phi_n(x) is x; evaluating it in 
                                % all the points gives us matrix "points"; 
                                % numb_pts-by-d     



 for j = 3:m_i_max              % For polynomial bases "phi_n", from n=3 to
                                % n=m_i_max, ...
    phi(:,:,j) = 2.* phi(:,:,2).*phi(:,:,j-1) - phi(:,:,j-2); 
                                % Use the recurrence formula to compute the 
                                % Chebyshev polynomial basis functions of 
                                % the degrees from 2 to m_i_max-1
 end 

% 1a. Rudik addition: compute non-mixed derivatives up to order 3
% Compute derivatives of unidimensional Chebyshev polynomials at these points
% d/dx T_n(x) = n*U_{n-1}(x)
if any(deriv > 0)
    
    deriv_index1 = find(deriv == 1);
    deriv_index2 = find(deriv == 2);
    deriv_index3 = find(deriv == 3);
    
    % Calculate polynomials of the second kind for derivative recurrence
    Un = ones(numb_pts,d,m_i_max);
    Un(:,:,2) = 2*points;
    
    for j = 3:m_i_max
        
        % U_{t}(x) = 2*x*U_{t-1}(x) - U_{t-2}(x)
        Un(:,:,j) = 2.*points.*Un(:,:,j-1) - Un(:,:,j-2);
        
    end
    
    if any(deriv == 1)
        % First derivatives and cross-partials
        
        grad = zeros(numb_pts,d,m_i_max);
        grad(:,:,2) = ones(size(points));

        for j = 3:m_i_max  

            % d/dx T_n(x) = n*U_{n-1}(x)
            % Since array indices start at 1, we require j-1 instead of j for
            % n, since indexing in grad and Un is relative, we just require
            % index j of grad and index j-1 of Un
            grad(:,:,j) = (j-1)*Un(:,:,j-1);          

        end

    end
    
    % Second derivatives
    % d^2Tn(x)/dx^2 = n*( (n+1)Tn(x) - Un(x)) / (x^2-1) ) if x == +/-1
    % d^2Tn(x)/dx^2 = (n^4 - n^2)/3 if x == 1
    % d^2Tn(x)/dx^2 = (-1)^n * (n^4 - n^2)/3 if x == -1
    if any(deriv ==2)
        
        % 0 for n = 0,1; 4 for n=2
        grad2 = zeros(numb_pts,d,m_i_max);
        grad2(:,:,3) = 4*ones(size(points));
        
        
        % For n = 3,...,m_i_max
        for j = 4:m_i_max  
            
            grad2(:,:,j) = (j-1)*( (j*phi(:,:,j) - Un(:,:,j))./(points.^2-1));
            if any(isnan(grad2(:,deriv_index2,j)))
                 if points(:,deriv_index2) == 1
                    
                    grad2(:,deriv_index2,j) = ((j-1)^4-(j-1)^2)/3;
                    
                elseif points(:,deriv_index2) == -1
                    
                    grad2(:,deriv_index2,j) = (-1)^(j-1)*(((j-1)^4-(j-1)^2)/3);
                    
                else
                    
                    error('Shouldnt get an infinity if the contracted point is not at 1.');
                    
                end
            end
        end
        
            
        
    end
    
    % Third derivatives
    % d^3Tn(x)/dx^3 = slightly more complicated
    if any(deriv == 3)
        
        % Need one more polynomial for third derivative
        phi(:,:,m_i_max+1) = 2.* phi(:,:,2).*phi(:,:,m_i_max) - phi(:,:,m_i_max-1); 
        
        % 0 for n = 0,1; 4 for n=2, 24 for n=3
        grad3 = zeros(numb_pts,d,m_i_max);
        grad3(:,:,4) = 24*ones(size(points));
        
        % For n = 4,...,m_i_max
        for j = 5:m_i_max
            
            grad3(:,:,j) = (j-1)*( ( (j-1+1)*(j-1).*Un(:,:,j-1) - ...
                ((j-1+1).*phi(:,:,j+1)-points.*Un(:,:,j))./(points.^2-1)).*(points.^2-1) - ...
                (2.*points.*( (j-1+1).*phi(:,:,j) - Un(:,:,j) ) ) )...
                ./(points.^2-1).^2;
            if any(isnan(grad3(:,deriv_index3,j)))
                if points(:,deriv_index3) == 1
                    
                    grad3(:,deriv_index3,j) = ((j-1)^4-(j-1)^2)*((j-1)^2-4)/15;
                    
                elseif points(:,deriv_index3) == -1
                    
                    grad3(:,deriv_index3,j) = (-1)^(3+j-1)*(((j-1)^4-(j-1)^2*((j-1)^2-4))/3);
                    
                else
                    
                    error('Shouldnt get an infinity if the contracted point is not at 1.');
                    
                end
            end
            
        end
        % Remove additional polynomial
        phi(:,:,m_i_max+1) = [];
        
    end
    
    % Replace unidimensional polynomials with derivatives
    if any(deriv == 1)
        phi(:,deriv_index1,:) = bsxfun(@rdivide,grad(:,deriv_index1,:),dist(deriv_index1));
    end
    if any(deriv == 2)
        phi(:,deriv_index2,:) = bsxfun(@rdivide,grad2(:,deriv_index2,:),dist(deriv_index2).^2);
    end
    if any(deriv == 3)
        phi(:,deriv_index3,:) = bsxfun(@rdivide,grad3(:,deriv_index3,:),dist(deriv_index3).^3);
    end

end
 
% 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
% required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
% and 3.4.2 for examples
% ----------------------------------------------------------------------


                                    % Compute the number of terms (i.e., multi-
numb_terms = size(Smol_elem,1);     % dimensional polynomial bases) in Smolyak 
                                    % polynomial     

Smol_bases = zeros(numb_pts,numb_terms);  % Initially, the matrix of multidimensional 
                                          % polynomial bases is empty
                                            
 
 for jt = 1:numb_terms         % For each term of Smolyak polynomial, ...
                                
     index_row = Smol_elem(jt,:);
     
     if any(deriv > 0)
         
         % Dimension of the polynomial
         poly_dim = numel(index_row(index_row > 1));
         
                                   % Identify the subindices of the unidimensional
                                   % basis function that constitute a jt-th multi-
                                   % dimensional basis function; this is a jt-th
                                   % row of matrix "Smol_elem"; 1-by-d
         if poly_dim == 0   
             
             product = zeros(numb_pts,1); % If there are no non-unity basis functions
                                          % the product is zero
                                         
         else
                                                   % If there is just one, check if it's the
             if any(deriv(find(index_row > 1)))    % dimension we're taking a derivative
                 product = ones(numb_pts,1);       % of. If it is, product may be non-zero.
             else 
                 product = zeros(numb_pts,1);
             end
             
         end
     
     
     else
         product = ones(numb_pts,1); 
     end
     
     
                               % Initialize a vector of products of unidimen-
                               % sional basis functions; numb_pts-by-1
     for jd = 1:d              % For each dimension (state variable), ...
         n = index_row(jd);    % A subindex of a unidimensional basis function 
                               % phi_n in a dimension jd is denoted n
         if n ~= 1;            % If the subindex of unidimensional basis 
                               % function is not equal to unity, ...
         product = product.*phi(:,jd,n); % USE PROD() HERE SOMEHOW?
                               % Compute the product of basis functions
                               
         % Otherwise, i.e., if n = 1, there is no need to compute the 
         % product of unidimensional basis functions, as it's equal to unity 
        end
     end

     Smol_bases(:,jt) = product; 
                               % Attach to the previously obtained matrix of 
                               % multidimensional basis functions a new
                               % product of unidimensional basis functions;
                               % e.g., for mu=1 and d=2, basis_bs is of
                               % size numb_pts-by-5
                               
 end
 
 

 
 