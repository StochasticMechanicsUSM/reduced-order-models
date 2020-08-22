%--------------------------------------------------------------------------
%Description: Function that computes the nw lowest frequencies and the 
%             corresponing mode shapes
%--------------------------------------------------------------------------
%Name: Unreduced_ModelAnalysis
%--------------------------------------------------------------------------
%Input:     - M     : Unreduced mass matrix                  (Double array)
%           - K     : Unreduced stiffness matrix             (Double array)
%           - Nw     : Number of mode shapes                (Integer array)
%--------------------------------------------------------------------------
%Output:    - w     : The Nw lowest frequencies              (Double array)
%           - PhiR  : The mode shapes corresponding to the Nw lowest
%                      frequencies                           (Double array)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [w,Phi] = Unreduced_ModelAnalysis(M,K,nw)
% Model dimension
NR = size(M,1);
if  nw <= NR
    % Eigenvalue problem
    options.v0 = ones(NR,1);
    [V,D,~]  = eigs(K,M,nw,'SM',options); 
    % Frequencies and mode shapes  
    [w,pos] = sort(sqrt(diag(D)));
    % Mode shapes at original physical coordinates
    Phi  = V(:,pos);  
else
    fprintf(['error: The number of frequencies nw must be less than '...
             num2str(NR) '\n'])
    w   = NaN;
    Phi = NaN;
end
return