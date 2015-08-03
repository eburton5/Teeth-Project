%% Goal of this script: calculate the cP distance between submeshes (size 400 pts) from genus Saimiri

%% Preparation
clearvars; 
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath([pwd '../Saimiri/off']));
% addpath(path,genpath([pwd '../MeshLabfiles/Saimiri/Boliviensis']));
% addpath(path,genpath([pwd '../MeshLabfiles/Saimiri/Sciureus']));
% addpath(path,genpath([pwd '../MeshLabfiles/Saimiri_Giseok']));

%%  Set Parameters for later
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 150;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

%% Prepare meshes
% Saimiri = cell(1,15);
% % Load Boliviensis species
% % MeshLabfiles/Saimiri/Boliviensis
% Saimiri{1} = Mesh('off','../Saimiri/off/aj02.off');
% Saimiri{2} = Mesh('off','../Saimiri/off/aj01.off'); 
% Saimiri{3} = Mesh('off','../Saimiri/off/aj04.off');
% Saimiri{4} = Mesh('off','../Saimiri/off/aj07.off');
% Saimiri{5} = Mesh('off','../Saimiri/off/aj09.off'); % replacement for aj03 (I cropped in Boyer lab and processed in Meshlab)
% % Saimiri{5} = Mesh('off','../Saimiri/off/aj03.off'); % this mesh has some significant artifacts 
% Saimiri{6} = Mesh('off','../Saimiri/off/as04.off');
% Saimiri{7} = Mesh('off','../Saimiri/off/aj06.off');
% % Load Sciureus species
% Saimiri{8} = Mesh('off','../Saimiri/off/ar22.off');
% Saimiri{9} = Mesh('off','../Saimiri/off/as01.off');
% Saimiri{10} = Mesh('off','../Saimiri/off/as03.off');
% Saimiri{11} = Mesh('off','../Saimiri/off/as02.off');
% Saimiri{12} = Mesh('off','../Saimiri/off/ar24.off');
% Saimiri{13} = Mesh('off','../Saimiri/off/ar23.off');
% Saimiri{14} = Mesh('off','../Saimiri/off/ar21.off');
% Saimiri{15} = Mesh('off','../Saimiri/off/aj05.off');
% 
% % Convert to .mat files
% for j = 1:15
%     [Saimiri{j}.Aux.Area,Saimiri{j}.Aux.Center] = Saimiri{j}.Centralize('ScaleArea');
%     Saimiri{j}.ComputeMidEdgeUniformization(options); % beginning of midedge uniformization addition
%     Saimiri{j}.Nf = Saimiri{j}.ComputeFaceNormals;
%     Saimiri{j}.Nv = Saimiri{j}.F2V'*Saimiri{j}.Nf';
%     Saimiri{j}.Nv = Saimiri{j}.Nv'*diag(1./sqrt(sum((Saimiri{j}.Nv').^2,1)));
%     Saimiri{j}.Aux.LB = Saimiri{j}.ComputeCotanLaplacian; 
%     %Saimiri{j} = Mesh(Saimiri{j}); % end of midedge uniformization addition
% end
% save('Saimiri_raw.mat', 'Saimiri'); %these can be used forever
% 
% % delete isolated vertices
% for j=1:15
%     ConfMax = Saimiri{j}.V(:,Saimiri{j}.Aux.ConfMaxInds);
%     GaussMax = Saimiri{j}.V(:,Saimiri{j}.Aux.GaussMaxInds);
%     GaussMin = Saimiri{j}.V(:,Saimiri{j}.Aux.GaussMinInds);
%     ADMax = Saimiri{j}.V(:,Saimiri{j}.Aux.ADMaxInds);
%     
%     dVInds = Saimiri{j}.DeleteIsolatedVertex();
%     if ~isempty(dVInds) %% recompute uniformization
%         [Saimiri{j}.Aux.Area,Saimiri{j}.Aux.Center] = Saimiri{j}.Centralize('ScaleArea');
%         [Saimiri{j}.BV, Saimiri{j}.BE] = Saimiri{j}.FindBoundaries();
%         ConfMax = (ConfMax-repmat(Saimiri{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Saimiri{j}.Aux.Area);
%         GaussMax = (GaussMax-repmat(Saimiri{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Saimiri{j}.Aux.Area);
%         GaussMin = (GaussMin-repmat(Saimiri{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Saimiri{j}.Aux.Area);
%         ADMax = (ADMax-repmat(Saimiri{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Saimiri{j}.Aux.Area);
%         
%         Saimiri{j}.ComputeMidEdgeUniformization(options);
%         Saimiri{j}.Nf = Saimiri{j}.ComputeFaceNormals;
%         Saimiri{j}.Nv = Saimiri{j}.F2V'*Saimiri{j}.Nf';
%         Saimiri{j}.Nv = Saimiri{j}.Nv'*diag(1./sqrt(sum((Saimiri{j}.Nv').^2,1)));
%         Saimiri{j}.Aux.LB = Saimiri{j}.ComputeCotanLaplacian;
%         
%         TREE = kdtree_build(Saimiri{j}.V');
%         Saimiri{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
%         Saimiri{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
%         Saimiri{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
%         Saimiri{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
%     end
% end
% 
% save('Saimiri_deliv.mat','Saimiri');

%% Using Giseok's meshes
% load Boliviensis species
% Saimiri{1} = load('../MeshLabfiles/Saimiri_Giseok/AMNH-M-38792_M904.mat'); Saimiri{1} = Saimiri{1}.G;
% Saimiri{2} = load('../MeshLabfiles/Saimiri_Giseok/AMNH-M-710003_M906.mat'); Saimiri{2} = Saimiri{2}.G;
% Saimiri{3} = load('../MeshLabfiles/Saimiri_Giseok/AMNH-M-76583_M908.mat'); Saimiri{3} = Saimiri{3}.G;
% Saimiri{4} = load('../MeshLabfiles/Saimiri_Giseok/AMNH-M-76586_M910.mat'); Saimiri{4} = Saimiri{4}.G;
% Saimiri{5} = load('../MeshLabfiles/Saimiri_Giseok/AMNH-M-98272_M912.mat'); Saimiri{5} = Saimiri{5}.G;
% Saimiri{6} = load('../MeshLabfiles/Saimiri_Giseok/USNM-364497_M914.mat'); Saimiri{6} = Saimiri{6}.G;
% Saimiri{7} = load('../MeshLabfiles/Saimiri_Giseok/USNM-396265_M916.mat'); Saimiri{7} = Saimiri{7}.G;
% % load Sciureus species
% Saimiri{8} = load('../MeshLabfiles/Saimiri_Giseok/USNM-518547_M918.mat'); Saimiri{8} = Saimiri{8}.G;
% Saimiri{9} = load('../MeshLabfiles/Saimiri_Giseok/USNM-545893_M920.mat'); Saimiri{9} = Saimiri{9}.G;
% Saimiri{10} = load('../MeshLabfiles/Saimiri_Giseok/USNM-546267_M922.mat'); Saimiri{10} = Saimiri{10}.G;

% load('Saimiri.mat');
% NumPts = 400;
% rng('shuffle');
% 
% % preallocations for speed
% i = zeros(1,15);
% DSai = cell(1,15);
% subV = cell(1,15);
% subF = cell(1,15);
% Sai = cell(1,15);
% 
% for j = 1:15 
%     i(j) = randi(Saimiri{j}.nV); % generate random starting point, will be different every time
%     % get subsampled mesh
%     subV{j} = Saimiri{j}.GeodesicFarthestPointSampling(NumPts,i(j));
%     [subV{j},subF{j}] = PerformGeodesicDelauneyTriangulation(Saimiri{j},subV{j},[]); % Create new mesh for subsample
%     Sai{j} = Mesh('VF',subV{j},subF{j});
%     % compute uniformization for subsampled mesh
%     [Sai{j}.Aux.Area, Sai{j}.Aux.Center] = Sai{j}.Centralize('ScaleArea');
%     [~,TriAreas] = Sai{j}.ComputeSurfaceArea;
%     Sai{j}.Aux.VertArea = (TriAreas'*Sai{j}.F2V)/3;
%     DSai{j} = Sai{j}.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
% %     Sai{j}.ComputeMidEdgeUniformization(options); % beginning of midedge uniformization addition
% %     Sai{j}.Nf = Sai{j}.ComputeFaceNormals;
% %     Sai{j}.Nv = Sai{j}.F2V'*Sai{j}.Nf';
% %     Sai{j}.Nv = Sai{j}.Nv'*diag(1./sqrt(sum((Sai{j}.Nv').^2,1)));
% %     Sai{j}.Aux.LB = Sai{j}.ComputeCotanLaplacian; 
% %     Saimiri{j} = Mesh(Sai{j}); % end of midedge uniformization addition
%     Sai{j}.Aux.UniformizationV = DSai{j}.V;
% %%%     Estimate conformal factor using area distortion
%     [~,AG] = Sai{j}.ComputeSurfaceArea;
%     [~,DG] = DSai{j}.ComputeSurfaceArea;
%     Sai{j}.Aux.Conf = (AG./DG)'*DSai{j}.F2V/3;
% %%%     set the G.Aux.Conf to 0 on the boundary
%     [Sai{j}.BV,Sai{j}.BE] = Sai{j}.FindBoundaries();
%     Sai{j}.Aux.Conf(Sai{j}.BV) = 0;
% 
% %%%     keep features on the original fine mesh for the downsampled mesh
%     TREE = kdtree_build(Sai{j}.V');
%     Sai{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, Saimiri{j}.V(:,Saimiri{j}.Aux.ConfMaxInds)');
%     %% similarly, find nearest neighbors of
%     %%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
%     Sai{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, Saimiri{j}.V(:,Saimiri{j}.Aux.ADMaxInds)');
%     Sai{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, Saimiri{j}.V(:,Saimiri{j}.Aux.GaussMaxInds)');
%     Sai{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, Saimiri{j}.V(:,Saimiri{j}.Aux.GaussMinInds)');
% 
%     %%% finding density points
%     minds = [Sai{j}.Aux.GaussMaxInds;Sai{j}.Aux.GaussMinInds;Sai{j}.Aux.ConfMaxInds];
%     minds = unique(minds);
%     Sai{j}.Aux.DensityPnts = Sai{j}.GeodesicFarthestPointSampling(Sai{j}.nV,minds);
% end
% save('Sai_400.mat','Sai'); 

%% delete isolated vertices
% for j=1:15
%     ConfMax = Sai{j}.V(:,Sai{j}.Aux.ConfMaxInds);
%     GaussMax = Sai{j}.V(:,Sai{j}.Aux.GaussMaxInds);
%     GaussMin = Sai{j}.V(:,Sai{j}.Aux.GaussMinInds);
%     ADMax = Sai{j}.V(:,Sai{j}.Aux.ADMaxInds);
%     
%     dVInds = Sai{j}.DeleteIsolatedVertex();
%     if ~isempty(dVInds) %% recompute uniformization
%         [Sai{j}.Aux.Area,Sai{j}.Aux.Center] = Sai{j}.Centralize('ScaleArea');
%         [Sai{j}.BV, Sai{j}.BE] = Sai{j}.FindBoundaries();
%         ConfMax = (ConfMax-repmat(Sai{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Sai{j}.Aux.Area);
%         GaussMax = (GaussMax-repmat(Sai{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Sai{j}.Aux.Area);
%         GaussMin = (GaussMin-repmat(Sai{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Sai{j}.Aux.Area);
%         ADMax = (ADMax-repmat(Sai{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Sai{j}.Aux.Area);
%         
%         Sai{j}.ComputeMidEdgeUniformization(options);
%         Sai{j}.Nf = Sai{j}.ComputeFaceNormals;
%         Sai{j}.Nv = Sai{j}.F2V'*Sai{j}.Nf';
%         Sai{j}.Nv = Sai{j}.Nv'*diag(1./sqrt(sum((Sai{j}.Nv').^2,1)));
%         Sai{j}.Aux.LB = Sai{j}.ComputeCotanLaplacian;
%         
%         TREE = kdtree_build(Sai{j}.V');
%         Sai{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
%         Sai{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
%         Sai{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
%         Sai{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
%     end
% end
% 
% save('Sai_1500.mat','Sai');
% 

%% load all subsamples from all 15 teeth
G = cell(1,150);
cPDist_Sai_all = zeros(150);

load('aj02_1000.mat'); G(1,1:10) = Sai;
load('aj02_1000Rslt.mat'); cPDist_Sai_all(1:10,1:10) = cPDist;
load('aj01_1000.mat'); G(1,11:20) = Sai;
load('aj01_1000Rslt.mat'); cPDist_Sai_all(11:20,11:20) = cPDist;
load('aj04_1000.mat'); G(1,21:30) = Sai;
load('aj04_1000Rslt.mat'); cPDist_Sai_all(21:30,21:30) = cPDist;
load('aj07_1000.mat'); G(1,31:40) = Sai;
load('aj07_1000Rslt.mat'); cPDist_Sai_all(31:40,31:40) = cPDist;
load('aj09_1000.mat'); G(1,41:50) = Sai;
load('aj09_1000Rslt.mat'); cPDist_Sai_all(41:50,41:50) = cPDist;
load('as04_1000.mat'); G(1,51:60) = Sai;
load('as04_1000Rslt.mat'); cPDist_Sai_all(51:60,51:60) = cPDist;
load('aj06_1000.mat'); G(1,61:70) = Sai;
load('aj06_1000Rslt.mat'); cPDist_Sai_all(61:70,61:70) = cPDist;
load('ar22_1000.mat'); G(1,71:80) = Sai;
load('ar22_1000Rslt.mat'); cPDist_Sai_all(71:80,71:80) = cPDist;
load('as01_1000.mat'); G(1,81:90) = Sai;
load('as01_1000Rslt.mat'); cPDist_Sai_all(81:90,81:90) = cPDist;
load('as03_1000.mat'); G(1,91:100) = Sai;
load('as03_1000Rslt.mat'); cPDist_Sai_all(91:100,91:100) = cPDist;
load('as02_1000.mat'); G(1,101:110) = Sai;
load('as02_1000Rslt.mat'); cPDist_Sai_all(101:110,101:110) = cPDist;
load('ar24_1000.mat'); G(1,111:120) = Sai;
load('ar24_1000Rslt.mat'); cPDist_Sai_all(111:120,111:120) = cPDist;
load('ar23_1000.mat'); G(1,121:130) = Sai;
load('ar23_1000Rslt.mat'); cPDist_Sai_all(121:130,121:130) = cPDist;
load('ar21_1000.mat'); G(1,131:140) = Sai;
load('ar21_1000Rslt.mat'); cPDist_Sai_all(131:140,131:140) = cPDist;
load('aj05_1000.mat'); G(1,141:150) = Sai;
load('aj05_1000Rslt.mat'); cPDist_Sai_all(141:150,141:150) = cPDist;
save('Sai_all_150_1000.mat','G');

% %% Create cP distance matrix
load('Saimiri.mat');
rslt = cell(150);
for i = 1:150
    for j = 1:150
        if cPDist_Sai_all(i,j) == 0
            try
                rslt{i,j} = ComputeContinuousProcrustes4(G{i},G{j},options); 
                cPDist_Sai_all(i,j) = rslt{i,j}.cPdist;
            catch
                rslt{i,j} = ComputeContinuousProcrustes4(G{j},G{i},options); 
                cPDist_Sai_all(i,j) = rslt{i,j}.cPdist;
            end
        end
    end
end

save('Sai_all_150_400_Rslt.mat','cPDist_Sai_all');

%% Compute distance from each subsampled mesh to each full mesh
CPrslt = cell(150,15);
CPDself = zeros(150,15); 
for i = 1:150
    for j = 1:15
        if CPDself(i,j) == 0
            CPrslt{i,j} = ComputeContinuousProcrustes4(G{i},Saimiri{j},options); %each row represents a single subsample, each column a single full mesh
            CPDself(i,j) = CPrslt{i,j}.cPdist;
        end
    end
end
save('Sai_all_150_1000_CPDself.mat','CPDself');