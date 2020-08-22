%--------------------------------------------------------------------------
%Description: Function that computes the nw lowest frequencies and the 
%             corresponing mode shapes
%--------------------------------------------------------------------------
%Name: RO_ModelAnalysis
%--------------------------------------------------------------------------
%Input:     - MR     : Reduced-order mass matrix             (Double array)
%           - KR     : Reduced-order stiffness matrix        (Double array)
%           - TR     : Reduced-order transformation matrix   (Double array)
%           - TG     : Transformation matrix                (Integer array)
%           - Nw     : Number of mode shapes                 (Integer array)
%--------------------------------------------------------------------------
%Output:    - wR     : The Nw lowest frequencies             (Double array)
%           - PhiR   : The mode shapes corresponding to the Nw lowest
%                      frequencies                           (Double array)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [wR,PhiR] = RO_ModelAnalysis(MR,KR,TR,TG,Nw)
% Reduced-order model dimension
NR = size(MR,1);
if  Nw <= NR
    % Eigenvalue problem
    options.v0 = ones(NR,1);
    [VR,DR,~]  = eigs(KR,MR,Nw,'SM',options); 
    % Frequencies and mode shapes  
    [wR,posR] = sort(sqrt(diag(DR)));
    % Mode shapes at original physical coordinates
    PhiR  = TG*TR*VR(:,posR);  
else
    fprintf(['error: The number of frequencies Nw must be less than '...
             num2str(NR) '\n'])
    wR   = NaN;
    PhiR = NaN;
end
return