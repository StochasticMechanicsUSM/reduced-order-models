%% Substructure definition and nominal reduced-order model
run('Substructure_coupling_Marco2p.m')
%% Model updating
% Nominal model parameters 
Theta0 = [1.0;
          1.0];
% Number of model parameters
n_theta = size(Theta0,1); 
% Mass model parameter function for each substructure     
g_th = @(Theta)[1;...
                1];
% Stiffness model parameter function for each substructure
h_th = @(Theta)Theta;

% Set of substructures that depends on the model parameter Thetaj j=1,2
S1 = [1,2];
S2 = [3,4];

S_th = {S1;
      S2};  
% New model parameters Thetak 
Thetak = [0.85;
          0.95];
%% Reduced-order model updating
fprintf('MODEL EVALUATED AT Thetak: \n \n')
% Substructure matrices and dominant-fixed interface modes updating
fprintf('Time of model updating (Normal modes):')
tic
[Components_Upd] = CompUpdating(Components,g_th(Thetak),h_th(Thetak),S_th,n_theta);
% Module IV-I
[MR_Upd,KR_Upd,TR_Upd,PerRed_Upd] = RO_Model(Components_Upd);
t1_R_Upd = toc;
fprintf([num2str(t1_R_Upd) '[s] \n'])
fprintf('Time of modal analisys (Normal modes):')
tic
% Modal analysis
[wR_Upd,PhiR_Upd]   = RO_ModelAnalysis(MR_Upd,KR_Upd,TR_Upd,TG,Nw);
t2_R_Upd = toc;
fprintf([num2str(t2_R_Upd) '[s] \n \n'])


               