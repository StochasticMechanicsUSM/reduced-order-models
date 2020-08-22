function [M,K] = UnreducedMatrices(GeneralQuantities,Joints,Mats,Sections,Elements)
% Initialization of mass and stiffness matrices
M = sparse(double(GeneralQuantities.NDoF),double(GeneralQuantities.NDoF));
K = sparse(double(GeneralQuantities.NDoF),double(GeneralQuantities.NDoF));
% Beam/Column elements
if GeneralQuantities.NBeamCol>0
    [Mbc,Kbc] = UnreducedBeamColMatrices(GeneralQuantities,Joints,Mats,Sections,Elements);
    M        = M + Mbc;
    K        = K + Kbc;
end
% Enforce symmetry
M = (M + M')/2; 
K = (K + K')/2; 
end

