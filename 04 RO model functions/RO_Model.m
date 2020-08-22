%--------------------------------------------------------------------------
%Description: Function that generates an reduced-order models based on
%             dominant fixed-interface normal modes
%--------------------------------------------------------------------------
%Name: Ro_Model
%--------------------------------------------------------------------------
%Input:     - Components                                  (Structure array)
%--------------------------------------------------------------------------
%Output:    - MR     : Reduced-order mass matrix             (Double array)
%           - KR     : Reduced-order stiffness matrix        (Double array)
%           - TR     : Reduced-order transformation matrix   (Double array)
%           - PerRed : Percentage of reduction                     (Double)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [MR,KR,TR,PerRed] = RO_Model(Components)
% Reduced-order transformation matrix
TR = [blkdiag(Components.Subs.Phi_id),...
      blkdiag(Components.Subs.Psi_ib)*Components.Tt;...
      sparse(Components.NI,Components.Nid),...
      speye(Components.NI,Components.NI)];
% Reduced-order mass matrix
MR = [speye(Components.Nid,Components.Nid),...
       blkdiag(Components.Subs.hat_Mib)*Components.Tt;...
       Components.Tt'*blkdiag(Components.Subs.hat_Mib)',...
       Components.Tt'*blkdiag(Components.Subs.hat_Mbb)*Components.Tt];
% Enforced symmetry
MR = (MR + MR')/2;
% Reduced-order stiffness matrix
KR = [blkdiag(Components.Subs.Lambda),...
       sparse(Components.Nid,Components.NI);...
       sparse(Components.NI,Components.Nid),...
       Components.Tt'*blkdiag(Components.Subs.hat_Kbb)*Components.Tt];
% Enforced symmetry
KR = (KR + KR')/2;
% Percentage of reduction
PerRed = 100*(Components.Ni-Components.Nid)/(Components.Ni+Components.NI);
return