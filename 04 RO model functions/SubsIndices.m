%--------------------------------------------------------------------------
%Description: Function that identifies and orders internal and boundary 
%             degrees of freedom
%--------------------------------------------------------------------------
%Name: SubsMatrices
%--------------------------------------------------------------------------
%Input:     - Joints            (Structure array)
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
function [Components] = SubsIndices(Joints,Elements,Components)
%% Identification of internal and boundary nodes of each substructure
for iSub = 1:Components.Ns
    % Auxiliary beam/column indices
    iBC_aux = Components.Subs(iSub).iBeamCol;
    % Substructure nodes
    Nod_aux = reshape(Elements.BeamCol.Nodes(iBC_aux,:),[],1);
    Nod_aux = sort(Nod_aux);
    % Non-repeated substructure nodes
    Components.Subs(iSub).Nod  = unique(Nod_aux);
    Components.Subs(iSub).Nnod = length(Components.Subs(iSub).Nod);
end
%% Internal and boundary nodes
for iSub = 1:Components.Ns
    % Identify substructures without current substructure(iSub)
    Aux_Ind = setdiff(1:Components.Ns,iSub)';
    % Non-repeated substructure nodes (Between substructures)
    Nod_aux = intersect(Components.Subs(iSub).Nod,...
                         vertcat(Components.Subs(Aux_Ind).Nod));
    % Boundary nodes
    Components.Subs(iSub).Nodb  = Nod_aux;
    Components.Subs(iSub).Nnodb = length(Components.Subs(iSub).Nodb);
    % Internal nodes
    Components.Subs(iSub).Nodi  = setdiff(Components.Subs(iSub).Nod,Nod_aux);
    Components.Subs(iSub).Nnodi = length(Components.Subs(iSub).Nodi);
    % Re-ordered Nodes
    Components.Subs(iSub).Nod = [Components.Subs(iSub).Nodi;...
                                  Components.Subs(iSub).Nodb];
end
%% Internal and boundary degrees of freedom at substructure level
for iSub = 1:Components.Ns
    % Internal degrees of freedom
    DoFsi = Joints.Connectivity(Components.Subs(iSub).Nodi,:);
    % Boundary degrees of freedom
    DoFsb = Joints.Connectivity(Components.Subs(iSub).Nodb,:);
    % Substructure degrees of freedom
    DoFs  = [DoFsi;...
             DoFsb];
    % Auxiliar vectors
    DoFs_aux  = reshape(DoFs',1,6*Components.Subs(iSub).Nnod)';
    DoFsi_aux = reshape(DoFsi',1,6*Components.Subs(iSub).Nnodi)';
    DoFsb_aux = reshape(DoFsb',1,6*Components.Subs(iSub).Nnodb)';
    % Number of internal and boundary degrees of freedom
    Components.Subs(iSub).NDoF  = sum(unique(DoFs_aux)  > 0);
    Components.Subs(iSub).NDoFi = sum(unique(DoFsi_aux) > 0);
    Components.Subs(iSub).NDoFb = sum(unique(DoFsb_aux) > 0);
    % Internal degrees of freedom at substructure level (local numbering)
    aux = find(DoFsi_aux > 0);
    C1  = sparse(aux,1:length(aux),1,length(DoFsi_aux),length(aux));
    [~,iA,iB] = unique(DoFsi_aux'*C1,'stable');
    C2 = sparse(1:length(iB),iB,1,length(iB),length(iA));
    Components.Subs(iSub).DoFsi = reshape(C1*C2*(1:length(iA))',6,...
                                          Components.Subs(iSub).Nnodi)';
    % Boundary degrees of freedom at substructure level (local numbering)
    aux = find(DoFsb_aux > 0);
    C1 = sparse(aux,1:length(aux),1,length(DoFsb_aux),length(aux));
    [~,iA,iB] = unique(DoFsb_aux'*C1,'stable');
    C2 = sparse(1:length(iB),iB,1,length(iB),length(iA));
    Components.Subs(iSub).DoFsb = reshape(C1*C2*(Components.Subs(iSub).NDoFi...
                                          +(1:length(iA)))',6,...
                                          Components.Subs(iSub).Nnodb)';
    % Degrees of freedom at local substructure level (local numbering)
    Components.Subs(iSub).DoFs = [Components.Subs(iSub).DoFsi;...
                                  Components.Subs(iSub).DoFsb];
end
% Save the number of internal and boundary degrees of freedom of all substructures
Components.Ni = sum(vertcat(Components.Subs.NDoFi));
Components.Nb = sum(vertcat(Components.Subs.NDoFb));
return