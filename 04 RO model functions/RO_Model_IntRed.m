%--------------------------------------------------------------------------
%Description: Function that generates an reduced-order models based on
%             dominant fixed-interface normal modes and interface reduction
%--------------------------------------------------------------------------
%Name: RO_Model_IntRed
%--------------------------------------------------------------------------
%Input:     - Components                                  (Structure array)
%--------------------------------------------------------------------------
%Output:    - MRI    : Reduced-order mass matrix             (Double array)
%           - KRI    : Reduced-order stiffness matrix        (Double array)
%           - TRI    : Reduced-order transformation matrix   (Double array)
%           - PerRedI: Percentage of reduction                     (Double)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [MRI,KRI,TRI,PerRedI] = RO_Model_IntRed(Components)
% Number of considered interface modes
NIR = size(Components.GammaI,2);
% Reduced-order transformation matrix
TRI = [blkdiag(Components.Subs.Phi_id),...
      blkdiag(Components.Subs.Psi_ib)*Components.Tt*Components.GammaI;...
      sparse(Components.NI,Components.Nid),...
      Components.GammaI];
% Reduced-order mass matrix
MRI = [speye(Components.Nid,Components.Nid),...
        blkdiag(Components.Subs.hat_Mib)*Components.Tt*Components.GammaI;...
        Components.GammaI'*Components.Tt'*blkdiag(Components.Subs.hat_Mib)',...
        speye(NIR,NIR)];
% Enforced symmetry
MRI = (MRI + MRI')/2;
% Reduced-order stiffness matrix
KRI = [blkdiag(Components.Subs.Lambda),...
      sparse(Components.Nid,NIR);...
      sparse(NIR,Components.Nid),...
      Components.LambdaI];
% Enforced symmetry
KRI = (KRI + KRI')/2;
% Percentage of reduction
PerRedI = 100*(Components.Ni + Components.NI-...
             (Components.Nid + NIR))/...
             (Components.Ni+Components.NI);
return