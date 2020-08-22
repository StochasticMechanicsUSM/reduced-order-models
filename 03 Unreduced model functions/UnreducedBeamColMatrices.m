function [Mef,Kef] = UnreducedBeamColMatrices(GeneralQuantities,Joints,Mats,Sections,Elements)
% Initialization of mass and stiffness matrices
Mef = sparse(double(GeneralQuantities.NDoF),double(GeneralQuantities.NDoF));
Kef = sparse(double(GeneralQuantities.NDoF),double(GeneralQuantities.NDoF));
% Number of beam/column elemets
NBC = GeneralQuantities.NBeamCol;
for iBC = 1:NBC
    % Transformation from local to global coordinates
    LocalAxes = BCLocalAxes(Joints.SpatialCoord,Elements.BeamCol.Nodes(iBC,:));
    % Transformation matrix
    aa   = sparse(blkdiag(LocalAxes,LocalAxes,LocalAxes,LocalAxes));
    % Element nodes
    inod = Elements.BeamCol.Nodes(iBC,1);
    jnod = Elements.BeamCol.Nodes(iBC,2);
    % Degrees of freedom
    gg = [Joints.Connectivity(inod,:),Joints.Connectivity(jnod,:)];  
    % Active degrees of freedom
    cc = sparse(find(gg>0),gg(gg>0),1,12,double(GeneralQuantities.NDoF));
    % Spatial coordinates
    nod = [Joints.SpatialCoord(Elements.BeamCol.Nodes(iBC,1),:);...
           Joints.SpatialCoord(Elements.BeamCol.Nodes(iBC,2),:)];
    % Local mass and stiffness matrices
    [Mf,Kf] = BCLocalMatrices(Mats,Sections,Elements,nod,iBC);
    % Global mass and stiffness matrices
    Mef = Mef + cc'*aa'*Mf*aa*cc;
    Kef = Kef + cc'*aa'*Kf*aa*cc;
end
return
    