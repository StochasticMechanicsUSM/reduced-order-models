function [Mef,Kef] = BeamColMatrices(Substructure,Joints,Mats,Sections,Elements)
% Initialization of mass and stiffness matrices
Mef = sparse(Substructure.NDoF,Substructure.NDoF);
Kef = sparse(Substructure.NDoF,Substructure.NDoF);
% Number of beam/column elements
NBC = length(Substructure.iBeamCol);
for jBC = 1:NBC
    iBC = Substructure.iBeamCol(jBC);
    % Transformation from local to global coordinates
    LocalAxes = BCLocalAxes(Joints.SpatialCoord,Elements.BeamCol.Nodes(iBC,:));
    % Transformation matrix
    aa   = sparse(blkdiag(LocalAxes,LocalAxes,LocalAxes,LocalAxes));     
    % Element nodes
    inod = Substructure.Nod == Elements.BeamCol.Nodes(iBC,1);
    jnod = Substructure.Nod == Elements.BeamCol.Nodes(iBC,2);
    % Degrees of freedom
    gg = [Substructure.DoFs(inod,:),Substructure.DoFs(jnod,:)];
    % Active degrees of freedom
    cc = sparse(find(gg>0),gg(gg>0),1,12,Substructure.NDoF);
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
    