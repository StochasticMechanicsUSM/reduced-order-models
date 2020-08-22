%--------------------------------------------------------------------------
%Description: Function that obtains dominant fixed-interface normal modes 
%             and interface constraint modes
%--------------------------------------------------------------------------
%Name: SubsModes
%--------------------------------------------------------------------------
%Input:     - Components                                  (Structure array)
%           - Nid : Vector of the number of dominant fixed-interface normal
%                   modes of each substructure             (Integer vector)
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
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [Components] = SubsModes(Components,Nid)
%% Substructure modes
% Total number of fixed-interface normal modes of all substructures
Components.Nid = sum(Nid);
for iSub = 1:Components.Ns
    if Nid(iSub) <= Components.Subs(iSub).NDoFi 
        % Number of fixed-interface normal modes per substructure
        Components.Subs(iSub).Nid = Nid(iSub);
        % Eigenvalues and dominant fixed-interface normal modes
        % Eigenproblem
        options.v0 = ones(Components.Subs(iSub).NDoFi,1);
        [V,D,~]    = eigs(Components.Subs(iSub).Kii,Components.Subs(iSub).Mii,...
                       Components.Subs(iSub).Nid,'SM',options);
        % Sort eigenvalues and the corresponding eigenvectors
        [D,pos] = sort(diag(D));
        V       = V(:,pos);
        % Normalization with respect to the internal mass matrix
        for i=1:Components.Subs(iSub).Nid
            V(:,i) = V(:,i)/sqrt(V(:,i)'*Components.Subs(iSub).Mii*V(:,i));
        end
        Components.Subs(iSub).Lambda = sparse(diag(D));
        Components.Subs(iSub).Phi_id = sparse(V);
        % Interface constraint modes
        Components.Subs(iSub).Psi_ib = -Components.Subs(iSub).Kii\...
                                        Components.Subs(iSub).Kib;
    else
        error(['Error: The number of dominant fixed-interface normal '... 
                 'modes of the substructure ' num2str(iSub) ' must be less '...
                 'than or equal to ' num2str(Components.Subs(iSub).NDoFi)])
    end
end
%% Auxiliarly substructure matrices
for iSub = 1:Components.Ns
   % Auxiliary mass matrices 
    Components.Subs(iSub).hat_Mib = Components.Subs(iSub).Phi_id'*...
                    Components.Subs(iSub).Mii*...
                    Components.Subs(iSub).Psi_ib + ...
                    Components.Subs(iSub).Phi_id'*Components.Subs(iSub).Mib;
    Components.Subs(iSub).hat_Mbb = Components.Subs(iSub).Psi_ib'*...
                    Components.Subs(iSub).Mii*...
                    Components.Subs(iSub).Psi_ib + ...
                    Components.Subs(iSub).Mib'*Components.Subs(iSub).Psi_ib +...
                    Components.Subs(iSub).Psi_ib'*Components.Subs(iSub).Mib +...
                    Components.Subs(iSub).Mbb;
   % Auxiliarly stiffness matrix
    Components.Subs(iSub).hat_Kbb = Components.Subs(iSub).Kib'*...
                    Components.Subs(iSub).Psi_ib +...
                    Components.Subs(iSub).Kbb; 
end
return