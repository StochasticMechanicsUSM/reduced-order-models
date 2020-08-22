function [Ms,Ks] = SubsMatricesGeneration(Substructure,Joints,Mats,Sections,Elements)
% Initialization of  mass and stiffness matrices
Ms = sparse(Substructure.NDoF,Substructure.NDoF);
Ks = sparse(Substructure.NDoF,Substructure.NDoF);
% Beam/Column elements
if ~isempty(Substructure.iBeamCol)
    [Mbc,Kbc] = BeamColMatrices(Substructure,Joints,Mats,Sections,Elements);
    Ms        = Ms + Mbc;
    Ks        = Ks + Kbc;
end
% Enforced symmetry
Ms = (Ms + Ms')/2;
Ks = (Ks + Ks')/2;
return

