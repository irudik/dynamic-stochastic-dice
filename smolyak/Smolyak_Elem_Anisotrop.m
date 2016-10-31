% Smolyak_Elem_Anisotrop.m is a routine that selects a subset of the subindices 
% of Smolyak elements corresponding to the given anisotropic case from a set of 
% subindices of the Smolyak isotropic elements; see "Smolyak method for solving  
% dynamic economic models: Lagrange interpolation, anisotropic grid and   
% adaptive domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and  
% Rafael Valero, (2014), Journal of Economic Dynamics and Control 44, 92–123 
% (henceforth, JMMV (2014)), Section 2.2.3 
%
% This version: Novenber 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Input:   "Smol_elem_iso"         is the matrix of subindices of unidi- 
%                                  mensional elements that constitute multi-
%                                  dimensional elements for the isotropic case
%          "vector_mus_dimensions" is the vector of the levels of
%                                  approximations in all dimensions  
% Output:  "Smol_elem_ani"         is the matrix of subindices of unidimen-
%                                  sional elements that constitute multidi-
%                                  mensional elements for the given anisotropic 
%                                  case
% 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


function  [Smol_elem_ani] = Smolyak_Elem_Anisotrop(Smol_elem_iso,vector_mus_dimensions)

 length_vmd = size(vector_mus_dimensions,2);% Compute the number of columns 
                                            % of the vector of mus across 
                                            % dimensions (is equal to the 
                                            % number of state variables)

 points_dimensions = zeros(1,length_vmd);   % This vector will tell how many
                                            % unidimensional elements in
                                            % each dimension we consider

 for i = 1:length_vmd
   aux = vector_mus_dimensions(1,i);
     if aux == 0                            % If the approximation level in
                                            % the i-th dimension is 0, ...
        points_dimensions(1,i) = 1;         % The number of unidimensional 
                                            % elements is 1
     else                                   % If the approximation level in
                                            % the i-th dimension is not 0,...
        points_dimensions(1,i) = 2^(aux)+1; % Compute the number of unidimensional
                                            % elements using the formula
     end
 end

 for i = 1:length_vmd
    aux1 = Smol_elem_iso(:,i) > points_dimensions(1,i);
    % If a subindex (i.e., number of elements) of the isotropic case is 
    % larger than that of the anisotropic case, set an auxiliary variable
    % "aux1" to 1
    Smol_elem_iso(aux1,: ) = [];  
    % Delete a vector of subindices for which aux=1
 end
Smol_elem_ani = Smol_elem_iso; 