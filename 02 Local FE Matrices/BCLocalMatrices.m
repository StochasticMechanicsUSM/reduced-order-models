function [Mf,Kf] = BCLocalMatrices(Mats,Sections,Elements,nod,iBC)
% Material and element properties
E         = Mats.E(Sections.BeamCol.iMat(Elements.BeamCol.iSecProp(iBC)));
U         = Mats.U(Sections.BeamCol.iMat(Elements.BeamCol.iSecProp(iBC)));
UMass     = Mats.UnitMass(Sections.BeamCol.iMat(Elements.BeamCol.iSecProp(iBC)));
AA        = Sections.BeamCol.SecProp.Area(Elements.BeamCol.iSecProp(iBC));
I33       = Sections.BeamCol.SecProp.I33(Elements.BeamCol.iSecProp(iBC));
I22       = Sections.BeamCol.SecProp.I22(Elements.BeamCol.iSecProp(iBC));
TorsConst = Sections.BeamCol.SecProp.TorsConst(Elements.BeamCol.iSecProp(iBC));
kap       = Sections.BeamCol.SecProp.AS2(Elements.BeamCol.iSecProp(iBC))/AA;
% Local mass and stiffness matrices
Mf = sparse(m_frame3D_mrHRZ1(UMass,AA,TorsConst,nod));
Kf = sparse(k_frame3D_ki(E,U,AA,I22,I33,TorsConst,kap,nod));
return
    