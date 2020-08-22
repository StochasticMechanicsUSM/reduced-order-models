%--------------------------------------------------------------------------
%Description: Function that generates transformation matrices that relate
%             boundary and interface degres of freedom and reorder the
%             physical degrees of freedom according to the original full 
%             finite element model
%--------------------------------------------------------------------------
%Name: TransMatrices
%--------------------------------------------------------------------------
%Input:     - Joints            (Structure array)
%           - Components        (Structure array)
%--------------------------------------------------------------------------
%Output:    - Components  : Structure Array that includes the substructure
%                           matrices                      (Structure array)
%           - TG          : Transformation matrix           (Integer array)
%--------------------------------------------------------------------------
%Notes:     - 
%--------------------------------------------------------------------------
%Assumptions: - Linear elastic models
%             - Substructuring setting was previously defined
%--------------------------------------------------------------------------
% Based on V. Araya's Work
% Written Nov 11 2019 by F. Mayorga franco.mayorga@gmail.com
%--------------------------------------------------------------------------
function [Components,TG] = TransMatrices(Joints,Components)
%% Transformation matrix
% All boundary nodes
Nodb = vertcat(Components.Subs.Nodb);
% Boundary degrees of freedom at local substructure level
DoFsb             = vertcat(Components.Subs.DoFsb);
DoFsb(DoFsb > 0)  = 1;
% Independent boundary nodes and the corresponding degrees of freedom    
[Nodb_i,iA,~]  = unique(Nodb,'stable');
DoFsb_i        = DoFsb(iA,:);
% Number of independent boundary degrees of freedom 
NDoFsb_i = sum(sum(DoFsb_i));
% Independent boundary degrees of freedom matrix (Interface DoFs)
DoFsb_i               = reshape(DoFsb_i',1,6*length(Nodb_i))';
DoFsb_i(DoFsb_i == 1) = 1:NDoFsb_i;
DoFsb_i               = reshape(DoFsb_i',6,length(Nodb_i))';
% Dependent boundary degrees of freedom matrix
for j = 1:size(Nodb,1)
    iaux = Nodb(j) == Nodb_i;
    DoFsb(j,:) = DoFsb_i(iaux,:);
end
DoFsb_aux = reshape(DoFsb',1,6*length(Nodb))';
DoFsb_d   = DoFsb_aux(DoFsb_aux > 0);
% Transformation matrix
Components.Tt = sparse(1:Components.Nb,DoFsb_d,1,Components.Nb,NDoFsb_i);
% Save the number of independet interface degrees of freedom
Components.NI = NDoFsb_i;
%% Transformation matrix TG
% Internal nodes
Nodi  = vertcat(Components.Subs.Nodi);
% Nodes of all substructures
Nod_s = [Nodi;Nodb_i];
% Degrees of freedom of all substructures
DoFs_s    = Joints.Connectivity(Nod_s,:);
DoFs_s_re = reshape(DoFs_s',1,6*size(Joints.Connectivity,1))';
% Transformation matrix
aux = find(DoFs_s_re > 0);
C1  = sparse((1:length(aux))',aux,1,length(aux),length(DoFs_s_re));
% Transformation matrix TG
[DoFs_s,iA,~] = unique(C1*DoFs_s_re,'stable');
TG            = sparse(DoFs_s,(1:length(iA))',1,length(DoFs_s),length(iA));
return