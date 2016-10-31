% Smolyak_Grid.m is a routine that constructs a multidimensional Smolyak  
% grid in the hypercube [-1,1]^d; see "Smolyak method for solving dynamic 
% economic models: Lagrange interpolation, anisotropic grid and adaptive 
% domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero,
% (2014), Journal of Economic Dynamics and Control 44, 92–123 (henceforth, 
% JMMV (2014)), Section 2.2.3 
%
% This version: Novenber 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Input:   "d"         is the number of dimensions (the number of state 
%                      variables)
%          "mu"        is the level of approximation (in the anisotropic
%                      case, this is maximum level of approximation)
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
% Output:  "Smol_grid" is the multidimensional Smolyak grid 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


function Smol_grid = Smolyak_Grid(d,mu,Smol_elem)

% 1. Compute the vector of extrema of Chebyshev polynomials corresponding 
% to the given level of Smolyak approximation mu
% -----------------------------------------------------------------------

% These points will be ordered as in Section 2.2.1 of JMMV(2014); e.g., for
% mu=1, the set of points is {0,-1,1}
         
points_1d = [];                    % Initially, the set of unidimensional 
                                   % points "points_1d" is empty; see JMMV  
                                   % (2014), Section 2.2.1
i_max = mu+1;                      % The maximum subindex of unidimensional
                                   % set A_i whose points are used to
                                   % construct Smolyak grid of the given mu; 
                                   %  e.g., for mu=1, we consider up to 
                                   % A_i_max={-1,1} where i_max=1+1=2
 for i = 1:i_max                   % A subindex of a unidimensional set of
                                   % points                                
     % Compute the number of elements, m(i),(using m(i)=2^(i-1)+1) in the  
     % i-th unidimensional set of points; see Section 2.2.1 in JMMV (2014)
     %---------------------------------------------------------------------
     m_i(i==1) = 1;                % If i=1, then m(i)=1
     m_i(i>1) =  2.^(i(i>1)-1) + 1;% If i>1, then m(i) = 2^(i-1)+1
                                   
     
     % Construct the extrema of Chebyshev polynomials used as unidimensional 
     % grid points in the Smolyak method
     %---------------------------------------------------------------------
     if (m_i==1)
        extrem_Cheb_1d = 0;
     else 
       
        j = (1:m_i)';                       
                                   % For j=1,...,m_i,...
        extrem_Cheb_1d = -cos(pi*(j-1)/(m_i-1));    
                                   % Chebyshev polynomials are defined in 
                                   % the interval [-1,1]
        extrem_Cheb_1d(abs(extrem_Cheb_1d)<1d-12) = 0;      
                                   % Round "extrem_Cheb_1d" to 0 if its  
                                   % absolute value is smaller than 1d-12
        extrem_Cheb_1d(1-extrem_Cheb_1d<1d-12) = 1;         
                                   % Round "extrem_Cheb_1d" to 1 if    
                                   % 1-extrem_Cheb_1d is smaller than 1d-12
        extrem_Cheb_1d(1+extrem_Cheb_1d<1d-12) = -1;        
                                   % Round "extrem_Cheb_1d" to -1 if   
                                   % 1+extrem_Cheb_1d is smaller than 1d-12
                                        
     end 
    
     points_1d = cat(1,points_1d, extrem_Cheb_1d);
                                   % Add to the previous set "points_1d" new 
                                   % points (given by extrema of unidimensional  
                                   % Chebyshev polynomials) as i increases
    
     points_1d = unique(points_1d,'rows','stable');
                                   % Choose the unrepeated points and order 
                                   % them as in Section 2.2.1 of JMMV (2014); 
                                   
    % !!! NOTICE: For the versions of MATLAB older than MATLAB R2012a,  
    % option 'stable' might not work
    
 end              

 
 % 2. Construct the matrix multidimensional Smolyak grid points for the   
 % required level of Smolyak approximation, mu; see JMMV (2014), Sections 2.2.3 
 % for examples
 % -------------------------------------------------------------------------
 Smol_grid = zeros(size(Smol_elem));
                               % Initialize the matrix of multidimensional
                               % Smolyak grid points
                               
 numb_points_md = size(Smol_grid,1);
                               % Compute the number of multidimensional 
                               % Smolyak grid points                                 
 
 for jp = 1:numb_points_md     % For each multidimensional grid point, ...
                                
     index_row = Smol_elem(jp,:);
                               % Identify the subindex of the unidimensional
                               % grid point jp; this is a jp-th row of matrix 
                               % "Smol_elem"; 1-by-d
    
     for jd = 1:d              % For each dimension (state variable), ...
         n = index_row(jd);    % A subindex of a unidimensional grid point 
                               % in a dimension jd is denoted n
        
         Smol_grid(jp,jd) = points_1d(n);
                               % Find the corresponding unidimensional grid
                               % point in the vector "points_1d"
     end 
 end
