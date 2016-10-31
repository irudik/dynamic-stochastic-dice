% Smolyak_Elem_Isotrop.m is a routine that constructs the subindices of the
% Smolyak elements (grid points and basis functions) for the isotropic case;
% see "Smolyak method for solving dynamic economic models: Lagrange interpo-
% lation, anisotropic grid and adaptive domain" by Kenneth L. Judd, Lilia Maliar, 
% Serguei Maliar and Rafael Valero, (2014), Journal of Economic Dynamics and 
% Control 44, 92–123 (henceforth, JMMV (2014)), Section 2.2.3
%
% This version: Novenber 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Input:   "d" is the number of dimensions (the number of state  variables) 
%          "mu" is the level of approximation
%  
% Output:  "Smolyak_elem_iso" is the vector of subindices of unidimensional 
%           elements (Smolyak grid points or polynomial basis function) that 
%           constitute a multidimensional element (Smolyak grid point or 
%           polynomial basis function) for the isotropic case
% 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

 function Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu)

% 1. Identify the indices of disjoint sets A's, i1,...,id, that satisfy the 
% Smolyak rule, d<=i1+i2+...+id<=|i|; see equation (1) in JMMV(2014)
% -------------------------------------------------------------------------
   
    Smol_rule = [];     % This will be a matrix of subindices i1,...,id of 
                        % disjoint sets A_i's such that their sum across all  
                        % dimensions, |i|=i1+i2+...+id,  satisfies  the 
                        % Smolyak rule, d<=|i|<=d+mu; initially, the matrix 
                        % is empty; at the intermediate steps j=0,...,mu, 
                        % the subindices satisfy d<=|i|=d+j
                       
    incr_Smol_rule = [];% The matrix of subindices i1,...,id of disjoint   
                        % sets A_i's, such that their sum across all dimensions,  
                        % |i|=i1+i2+...+id, is exactly equal to d+j, i.e., 
                        % |i|=d+j; initially, the matrix is empty; when j 
                        % increases, this matrix is concatinated to matrix  
                        % "Smol_rule" obtained for the previous value of j,
                        % i.e., j-1
   for j = 0:mu  
       
     prev_incr = incr_Smol_rule; 
                       % Ffor the previous value of j, call "incr_Smol_rule"  
                       % as "prev_incr"
     
     % Identify new subindices of unidimensional sets A's i1,i2,...,id that 
     % jointly satisfy the Smolyak rule as the sum of subindices increases 
     % from j-1 to j; see JMMV (2014), Section 2.2 
     %---------------------------------------------------------------------                
     if (j == 0)
       incr_Smol_rule = ones(1,d); % The matrix of subindices for j=0 is a 
                                   % 1-by-d vector of ones, (1,...,1)
     
     else  
       m = size(prev_incr,1);
                                   % Compute the number of columns in the
                                   % previous matrix of subindices that
                                   % satisfy the Smolyak rule "prev_incr"
       incr_Smol_rule = [];         % Initially, the matrix of subindices is 
                                   % empty   
       aux = zeros(m,d);           % Allocate memory to an initial auxiliary
                                   % matrix that will be added to
                                   % "prev_incr" 
        for id = 1:d
            aux_new = aux;         % New auxiliary matrix is equal to the old 
                                   % one
            aux_new(:,id) = 1;      % For a dimension i, set elements of this
                                   % new auxiliary matrix to 1 instead of 0
                                   
            augmented = prev_incr + aux_new; 
                                   % Increase the subinices of
                                   % "prevoius_incr" by 1 in dimension id
            
            incr_Smol_rule = cat(1,incr_Smol_rule, augmented);                       
                                   % Concatenate "incr_Smol_rule" and
                                   % "augmented": the first row of "augmented"
                                   % goes after the last row of
                                   % "incr_Smol_rule" 
        end
       
        
     end
         incr_Smol_rule = unique(incr_Smol_rule,'rows');
                                   % Eliminate the repeated indices in the
                                   % matrix "incr_Smol_rule"   
                                   
     % Concatenate the matrix of newly constructed indices to the previously 
     % obtained matrix of indices satisfying the Smolyak rule
     %---------------------------------------------------------------------
     Smol_rule = cat(1,Smol_rule, incr_Smol_rule); % E.g., for mu=1 and d=2,
                                                   % Smol_rule=[1 1; 1 2; 2 1]
                       
   end
  
   n_comb = size(Smol_rule,1);
              % The total number of combinations of indices, i1,...,id, 
              % of the unidimensional disjoint sets A's that satisfy the 
              % Smolyak rule; e.g., for mu=1 and d=2, n_comb=3 such as (1)
              % i1=1 and i2=1, (2) i1=1 and i2=2, (3) i1=2 and i2=1
              

% 2. Construct the multidimensional indices of elements as a Cartesian product
% of unidimensional indices
% -------------------------------------------------------------------------

Smolyak_elem_iso = []; 
                % The matrix of multidimensional indices of elements (points)
                %  belonging to the disjoint subsets A's that satisfy the 
                % Smolyak rule (every single point is indexed); initially, 
                % the matrix is empty
              
for i = 1:n_comb            % For each combination of subindices i1,...,id 
                            % that satisfies the Smolyak rule
    
   
    incr_indices = [];      % This is a matrix of multidimensional indices
                            % of unidimensional grid points that will be
                            % added to the final multidimensional matrix of 
                            % indices "Smolyak_elem_iso", as n_comb increases;
                            % initially, the matrix is empty

    one_comb = Smol_rule(i,:); 
                            % Consider the i-th row of "Smol_rule" (i.e.,
                            % one particular combination of the subindices 
                            % of the disjoint sets satisfying the Smolyak
                            % rule); e.g., for mu=1 and d=2, Smol_rule=
                            % [1 1; 1 2; 2 1] and one_comb for i=1 is [1 1];
                            % 1-by-d
    for jd = 1:d            % For each dimension jd, ...
        
      prev_indices = incr_indices;   
            
          % Compute the indices of elements (points) of the unidimensional 
          % set
          %----------------------------------------------------------------
            if one_comb(jd) == 1                        
                           % Take a jd-th element of the row vector "one_comb" 
                           % that corresponds to dimension jd (this is the 
                           % subindex of the unidimensional disjoint set A_i 
                           % from which this element comes from; if an element
                           % (point) is from the disjoint set
                           % A_1,...              
               indices_elem_jd = 1;                        
                           % A_1 contains one element; this element is indexed "1"
   
            elseif one_comb(jd) == 2                    
                           % If an element (point) is from the disjoint set
                           % A_2,...
               indices_elem_jd = [2;2^(one_comb(jd)-1)+1];            
                           % A_2 contains two elements; these elements are
                           % indexed "2" and "3", so that indices_elem_jd=[2;3]
            else
               indices_elem_jd = (2^(one_comb(jd)-2)+2 : 2^(one_comb(jd)-1)+1)'; 
                           % The subsequent disjoint sets contain the elements 
                           % from m(one_comb(jd)-1)+1=2^(one_comb(jd)-2)+2 to 
                           % m(one_comb(jd))=2^(one_comb(jd)-1)+1; e.g., the 
                           % disjoint set A_4 contains the elements indexed 
                           % 6, 7, 8 and 9, so that indices_elem_jd=[6,7,8,9]
            end                      
                           
           % Create a Cartesian product of two sets, "prev_indices" and 
           % "indices_elem_jd"
           %-----------------------------------------------------------
           a = prev_indices;   % Call one set "a" 
           b = indices_elem_jd;% Call the other set "b"
           if (isempty(b))     % If "b" is empty, ...
               z = a;          % The Cartesian product "z" is given by "a"
           elseif (isempty(a)) % If "a" is empty, ...
               z = b;          % The Cartesian product "z" is given by "b"   
           else                % For the case of non-empty sets, ...
               
               z = [];         % Initially, the Cartesian product is empty
               a_columns = size(a,1);  
               b_columns = size(b,1);
               for k = 1:b_columns;
                  z = [z;[a ones(a_columns,1)*b(k,:)]]; 
               end
           end
           
           incr_indices = z;   % 
    end
    
    Smolyak_elem_iso = cat(1,Smolyak_elem_iso,incr_indices); 
                           % Construct the matrix of multidimensional indices 
                           % of elements by concatenating "incr_indices" to 
                           % the previous "Smolyak_elem_iso" obtained for the 
                           % previous combination "n_comb"
end   

