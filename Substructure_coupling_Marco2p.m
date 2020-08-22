clear, clc, close all
addpath(genpath(pwd));
path(path)
% load model input information
load('marco_2p.mat')
% Reduced-order model generation
tic
% Module I
[Components] = SubsIndices(Joints,Elements,Components);
% Module II
[Components] = SubsMatrices(Joints,Mats,Sections,Elements,Components);
% Module III
[Components,TG] = TransMatrices(Joints,Components);
% Module IV
Nid = [3;3;1;1];
[Components] = SubsModes(Components,Nid);
t0_R = toc;
fprintf('Time of generation of RO model (Normal modes):')
tic
% Module V-I
[MR,KR,TR,PerRed] = RO_Model(Components);
t1_R = toc;
fprintf([num2str(t0_R + t1_R) '[s] \n'])
fprintf('Time of modal analisys of RO model (Normal modes):')
tic
% Modal analyses
Nw = 5;
[wR,PhiR]   = RO_ModelAnalysis(MR,KR,TR,TG,Nw);
t2_R = toc;
fprintf([num2str(t2_R) '[s] \n \n'])
fprintf('Time of generation of RO model (Interface reduction):')
tic
% Module V-II
NIR = 7;
[Components] = InterfaceModes(Components,NIR);
[MRI,KRI,TRI,PerRedI] = RO_Model_IntRed(Components);
t1_RI = toc;
fprintf([num2str(t0_R + t1_RI) '[s] \n'])
fprintf('Time of modal analisys of RO model (Interface reduction):')
tic
% Modal analyses
Nw = 5;
[wRI,PhiRI] = RO_ModelAnalysis(MRI,KRI,TRI,TG,Nw);
t2_RI = toc;
fprintf([num2str(t2_RI) '[s] \n \n'])
%% Unreduced model   
fprintf('Time of generation of unreduced model:')
tic
% Mass and stiffness unreduced matrices based on substructuring
[M,K] = UnreducedMatrices(GeneralQuantities,Joints,Mats,Sections,Elements);
t1 = toc;
fprintf([num2str(t1) '[s] \n'])
fprintf('Time of modal analisys of Unreduced model:')
tic
% Modal analyses
[w,Phi] = Unreduced_ModelAnalysis(M,K,Nw);  
t2 = toc;
fprintf([num2str(t2) '[s] \n\n'])