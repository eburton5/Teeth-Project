%% purpose of this script is to identify teeth that are likely cracked. Problems with cracked teeth arise when
% ConfMax points are placed on both sides of a crack, even though their is
% truly only one feature. This script will delete one of the two.

%% preparation
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath([pwd '/PNAS/meshes']));

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

%% load mesh and prepare variables for use
load('x23.mat'); %now have variable Gs
NumConfMax = length(Gs{1}.Aux.ConfMaxInds);
% ConfMaxDiff = cell(10,NumConfMax-1);
% ConfMaxDist = zeros(10,NumConfMax);
f = 0;

%% find distances between ConfMaxInds (this only found it for adjacent ConfMaxInds in list)
% for i = 1:10
%     for j = 2:NumConfMax
%         ConfMaxDiff{i,j-1} = Gs{i}.V(:,Gs{i}.Aux.ConfMaxInds(j)) - Gs{i}.V(:,Gs{i}.Aux.ConfMaxInds(j-1));
%     end
% end
% for i = 1:10
%     for j = 2:NumConfMax
%         ConfMaxDist(i,j-1) = sqrt(sum((Gs{i}.V(:,Gs{i}.Aux.ConfMaxInds(j)) - Gs{i}.V(:,Gs{i}.Aux.ConfMaxInds(j-1))).^2));
%     end
% end   

%% find distances between all pairs of ConfMaxInds; check for distance less than 0.1
for j = 1:NumConfMax
    for k = 1:NumConfMax
        ConfMaxDist = sqrt(sum((G.V(:,G.Aux.ConfMaxInds(k)) - G.V(:,G.Aux.ConfMaxInds(j))).^2));
        if ConfMaxDist == 0 
            continue
        elseif ConfMaxDist < 0.1
            disp(ConfMaxDist);
            f = j;
            break
        else
            continue
        end
    end
end
    
%% delete 1 confmax point if Euclidean distance between any two confmax points is less than 0.1
% tooclose = 0;
% 
% for i = 1:10
%     for j = 2:NumConfMax
%         if ConfMaxDist(i,j-1) < 0.1
% %         if min(abs(ConfMaxDiff{i,j-1})) < 0.015
%             tooclose = tooclose + 1;
%             f = j;
%         end 
%     end
% end

% if tooclose > 0
%     for k = 1:10
%        Gs{k}.Aux.ConfMaxInds(f) = [];
%     end
%     save('x23_fixed_for_crack.mat','Gs');
if f ~= 0
    G.Aux.ConfMaxInds(f) = [];
    save('x23_400_fixed_for_crack.mat','Gs');
    %% compute cPdist matrix
    cPDist = zeros(10);
    rslt = cell(10);
    for i = 1:10
        for j = 1:10
            try
                rslt{i,j} = ComputeContinuousProcrustes4(Gs{i},Gs{j},options); 
                cPDist(i,j) = rslt{i,j}.cPdist;
            catch
                rslt{i,j} = ComputeContinuousProcrustes4(Gs{j},Gs{i},options); 
                cPDist(i,j) = rslt{i,j}.cPdist;
            end
        end
    end
    save('x23_400_fixed_for_crack_Rslt.mat','cPDist');
end

