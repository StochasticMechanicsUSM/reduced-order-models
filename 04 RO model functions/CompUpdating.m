%--------------------------------------------------------------------------
%Description: Function that update modal quantities of each substructure
%             in terms of function g(theta) and h(theta)
%--------------------------------------------------------------------------
%Name: SubsModes
%--------------------------------------------------------------------------
%Input:     - Components                                  (Structure array)
%           - g : Vector of linear or nonlinear functions that modify the 
%                 substructure mass matrices in terms of the vector of
%                 model  parameters theta                   (Double vector)
%           - h : Vector of linear or nonlinear functions that modify the 
%                 substructure stiffness matrices in terms of the vector of 
%                 model parameters theta                    (Double vector)
%           - Sj: Set of substructures that depends on the model parameter 
%                 Thetaj                                       (Cell array)
%           - n_theta: Number of model parameters                 (Integer)
%--------------------------------------------------------------------------
%Output:    - Components : Structure Array that includes the fixed-
%                          interface normal modes and interface constraint 
%                          modes                          (Structure array)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga carlos.mayorga@usm.cl
%--------------------------------------------------------------------------
function [Components] = CompUpdating(Components,g_th,h_th,Sj,n_theta)
% Loop through the model parameters
for j = 1:n_theta    
    for Subj = Sj{j}
        %% Substructure matrices updating
        % Mass substructure matrices 
        if g_th(j)~=1 
            Components.Subs(Subj).Mii = Components.Subs(Subj).Mii*g_th(j);
            Components.Subs(Subj).Mib = Components.Subs(Subj).Mib*g_th(j);
            Components.Subs(Subj).Mbb = Components.Subs(Subj).Mbb*g_th(j);
        end
        % Stiffness substructure matrices
        if h_th(j)~=1 
            Components.Subs(Subj).Kii = Components.Subs(Subj).Kii*h_th(j);
            Components.Subs(Subj).Kib = Components.Subs(Subj).Kib*h_th(j);
            Components.Subs(Subj).Kbb = Components.Subs(Subj).Kbb*h_th(j);
        end
        %% Substructure modes updating
        % Dominant fixed-interface normal modes updating
        if g_th(j)~=1 
            Components.Subs(Subj).Phi_id = Components.Subs(Subj).Phi_id/...
                                           sqrt(g_th(j));
        end
        % Eigenvalues updating
        if g_th(j)~=1||h_th(j)~=1 
            Components.Subs(Subj).Lambda = Components.Subs(Subj).Lambda*...
                                           (h_th(j)/g_th(j));
        end
        %% Auxiliary substructure matrices updating
        % Auxiliary mass matrices 
        if g_th(j)~=1 
            Components.Subs(Subj).hat_Mib = Components.Subs(Subj).hat_Mib*...
                                            sqrt(g_th(j));
            Components.Subs(Subj).hat_Mbb = Components.Subs(Subj).hat_Mbb*...
                                            g_th(j);
        end
        % Auxiliary stiffness matrix
        if h_th(j)~=1 
            Components.Subs(Subj).hat_Kbb = Components.Subs(Subj).hat_Kbb*...
                                            h_th(j);
        end
    end
end
return