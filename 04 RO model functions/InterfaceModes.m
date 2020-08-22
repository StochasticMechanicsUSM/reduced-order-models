%--------------------------------------------------------------------------
%Description: Function that obtains the NIR kept interface modes and the 
%             corresponding eigen values
%--------------------------------------------------------------------------
%Name: RO_Model_IntRed
%--------------------------------------------------------------------------
%Input:     - Components                                  (Structure array)
%           - NI : Number of considered interface modes     (Integer array)
%--------------------------------------------------------------------------
%Output:    - Components : Structure Array that includes the interface 
%                          modes and their eigenvalues             (Double)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [Components] = InterfaceModes(Components,NIR)
% Interface reduction
if  NIR <= Components.NI
% With interface reduction
    % Interface mass and stiffness matrices
    MI = Components.Tt'*blkdiag(Components.Subs.hat_Mbb)*Components.Tt;
    KI = Components.Tt'*blkdiag(Components.Subs.hat_Kbb)*Components.Tt;
    % Interface frequencies and mode shapes
    options.v0 = ones(Components.NI,1);
    [VI,DI,~]  = eigs(KI,MI,NIR,'SM',options);
    % Sort eigenvalues and the corresponding eigenvectors
    [DI,pos_I]  = sort(diag(DI));
    VI          = VI(:,pos_I);
    % Normalization with respect to the interface mass matrix
    for i=1:NIR
        VI(:,i) = VI(:,i)/sqrt(VI(:,i)'*MI*VI(:,i));
    end
    % Interface modes and eigenvalues
    Components.GammaI  = sparse(VI);
    Components.LambdaI = sparse(diag(DI));
else
    error(['Error: The number of considered interface modes must be less than '...
           'or equal to ' num2str(Components.NI)])
end
return