%--------------------------------------------------------------------------
%Description: Function that generates substructure matrices
%--------------------------------------------------------------------------
%Name: SubsMatrices
%--------------------------------------------------------------------------
%Input:     - Joints            (Structure array)
%           - Mats              (Structure array)
%           - Sections          (Structure array)
%           - Elements          (Structure array)
%           - Components        (Structure array)
%--------------------------------------------------------------------------
%Output:    - Components  : Structure Array that includes the substructure
%                           matrices                      (Structure array)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [Components] = SubsMatrices(Joints,Mats,Sections,Elements,Components)
%% Generation of mass and stiffness matrices at substructure level
for iSub = 1:Components.Ns
    % Mass and stiffness matrices at substructure level
    [Components.Subs(iSub).Ms,Components.Subs(iSub).Ks] = ...
    SubsMatricesGeneration(Components.Subs(iSub),Joints,Mats,Sections,Elements);
    % Partition in terms of internal and boundaty degrees of freedom
    % Internal and boundary indices
    ind_i = 1:Components.Subs(iSub).NDoFi;
    ind_b = (1 + Components.Subs(iSub).NDoFi):Components.Subs(iSub).NDoF;
    % Internal and boundary mass matrices at substructure level
    Components.Subs(iSub).Mii = Components.Subs(iSub).Ms(ind_i,ind_i);
    Components.Subs(iSub).Mib = Components.Subs(iSub).Ms(ind_i,ind_b);
    Components.Subs(iSub).Mbb = Components.Subs(iSub).Ms(ind_b,ind_b);
    % Internal and boundary stiffness matrices at substructure level
    Components.Subs(iSub).Kii = Components.Subs(iSub).Ks(ind_i,ind_i);
    Components.Subs(iSub).Kib = Components.Subs(iSub).Ks(ind_i,ind_b);
    Components.Subs(iSub).Kbb = Components.Subs(iSub).Ks(ind_b,ind_b);
end
return