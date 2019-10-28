% Aggregate results from datagen_drive_alfa41ZPG into filtered snapshot
% lists associated with vortex core positions (y_over_delta and z_over_S)

%% Preparation
clear all; close all; clc; 

% Add Paths
addpath closure_relations/
addpath definitions/
addpath intersections/
addpath datagen/


%% Load snapshots from BLT ZPG straight case (Baldacchinno)
B001 = load('datagen/snapshotmats/B001.mat');

%% Load snapshots from Wendt 2001 (cases of table 5 ) (unusual flow orientation)
W001 = load('datagen/snapshotmats/W001.mat', 'snapshots');
W002 = load('datagen/snapshotmats/W002.mat', 'snapshots');
W003 = load('datagen/snapshotmats/W003.mat', 'snapshots');
W004 = load('datagen/snapshotmats/W004.mat', 'snapshots');

%% Load snapshots from Wendt 2001 (cases of Table 6 ) (  usual flow orientation)
W005 = load('datagen/snapshotmats/W005.mat', 'snapshots');
W006 = load('datagen/snapshotmats/W006.mat', 'snapshots');
W007 = load('datagen/snapshotmats/W007.mat', 'snapshots');
W008 = load('datagen/snapshotmats/W008.mat', 'snapshots');

%% Load snapshots from Logdberg 2009 (JFM paper 1st case)
L001 = load('datagen/snapshotmats/L001.mat', 'snapshots');
L002 = load('datagen/snapshotmats/L002.mat', 'snapshots');
L003 = load('datagen/snapshotmats/L003.mat', 'snapshots');


%% Now aggregate
% Create condensed snapshots object
CSN = condensed_snapshots();
% Add snapshots source by source!
CSN.add_snapshots(B001.snapshots)
% Table 5 of Wendt 2001 (unusual flow orientation)
CSN.add_snapshots(W001.snapshots)
CSN.add_snapshots(W002.snapshots)
CSN.add_snapshots(W003.snapshots)
CSN.add_snapshots(W004.snapshots)
% Table 6 of Wendt 2001 (  usual flow orientation)
CSN.add_snapshots(W005.snapshots)
CSN.add_snapshots(W006.snapshots)
CSN.add_snapshots(W007.snapshots)
CSN.add_snapshots(W008.snapshots)
% Logdberg 2009 (cases of JFM paper)
CSN.add_snapshots(L001.snapshots)
CSN.add_snapshots(L002.snapshots)
CSN.add_snapshots(L003.snapshots)

save('datagen/condensed_snapshots.mat', 'CSN')


%% Now sample the data

% % Filter out initial state (all 0's)
N_total_snapshots       = size(CSN.u_tilde_over_ue, 3);
% % Make selection of X snapshots
N_sampled_snapshots     = 200; %2023;
% Allocate X matrix
X = zeros(size(CSN.u_tilde_over_ue,1) * size(CSN.u_tilde_over_ue,2),N_sampled_snapshots); 
% Fill first column of X matrix
u_tilde_over_ue_sampled = CSN.u_tilde_over_ue(:,:,1); X(:,1) = u_tilde_over_ue_sampled(:);
% Fill next columns of X matrix
for n_sampled_snapshot = 2:N_sampled_snapshots
    n_ref_step = round((n_sampled_snapshot-1) * (N_total_snapshots / (N_sampled_snapshots-1)));
    u_tilde_over_ue_sampled = CSN.u_tilde_over_ue(:,:,n_ref_step);
    X(:,n_sampled_snapshot) = u_tilde_over_ue_sampled(:);
end

% Make SVD decomposition (economy size)
[Umat,Smat,Vmat] = svd(X, 0);


% Extract contribution of each mode to L2 norm from S matrix (note that det(V) = 1 and the columns of U are orthonormal)
lambda_S = zeros(1,N_sampled_snapshots);
for n_sampled_snapshot = 1:N_sampled_snapshots
    lambda_S(n_sampled_snapshot) = Smat(n_sampled_snapshot, n_sampled_snapshot);
end
% Now normalize
lambda_S_cumsum = cumsum(lambda_S) / sum(lambda_S);
% And plot
figure(10); 
subplot(224)
% plot(1:N_condensed_snapshot, lambda_S_cumsum); 
stairs(0:N_sampled_snapshots, [0 , lambda_S_cumsum]); 
hstairs = stairs(0:N_sampled_snapshots, [0 , lambda_S_cumsum], 'x-'); hold on;
plot(1, lambda_S_cumsum(1), 'o-');
plot(2, lambda_S_cumsum(2), 'o-');
plot(3, lambda_S_cumsum(3), 'o-');
grid on; axis([0 20 0 1]);
title('Normalized mode eigenvalues');
xlabel('N_{mode}'); ylabel('Eigenvalues of S matrix');
legend('Cumulative Sum (\rightarrow 1)', ['Scaled \lambda_1 = ' , num2str(lambda_S_cumsum(1))], ['Scaled \lambda_2 = ' , num2str(lambda_S_cumsum(2)-lambda_S_cumsum(1))], ['Scaled \lambda_3 = ' , num2str(lambda_S_cumsum(3)-lambda_S_cumsum(2))], 'Location', 'South')

% Extract first three (four) modes, and reshape them back into matrices
phi1_vec    = Umat(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = Umat(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = Umat(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = Umat(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot them
figure(101)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')

subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')

subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')

subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')


%% Now play with Laurens van der Maaten Dimensionality Reduction Toolbox
addpath drtoolbox/
addpath drtoolbox/techniques/

% Use PCA from DR Toolbox
[mappedX_PCA, mapping_PCA] = compute_mapping(X, 'PCA', 4);

% Try it out now
phi1_vec    = mappedX_PCA(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = mappedX_PCA(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = mappedX_PCA(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = mappedX_PCA(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot it!
figure(201)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')

subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')

subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')

subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')


%% Play with options : Autoencoder
techniques = {'PCA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 'Sammon', 'Isomap', ...
        'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', ...
        'FastMVU', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', ...
        'LLTSA', 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA'};
    
% Use Autoencoder from DR Toolbox
[mappedX_AE, mapping_AE] = compute_mapping(X', 'Autoencoder', 4);

% Try it out now
phi1_vec    = mappedX_AE(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = mappedX_AE(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = mappedX_AE(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = mappedX_AE(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot it!
figure(301)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')


subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')


subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')


subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')



%% Play with options : Kernel PCA
techniques = {'PCA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 'Sammon', 'Isomap', ...
        'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', ...
        'FastMVU', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', ...
        'LLTSA', 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA'};
    
% Use Autoencoder from DR Toolbox
[mappedX_KPCA, mapping_KPCA] = compute_mapping(X, 'KernelPCA', 4);

% Try it out now
phi1_vec    = mappedX_KPCA(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = mappedX_KPCA(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = mappedX_KPCA(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = mappedX_KPCA(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot it!
figure(401)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')


subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')


subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')


subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')

%% Play with options : SNE
techniques = {'PCA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 'Sammon', 'Isomap', ...
        'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', ...
        'FastMVU', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', ...
        'LLTSA', 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA'};

% Use Autoencoder from DR Toolbox
[mappedX_SNE, mapping_SNE] = compute_mapping(X, 'SNE', 4);

% Try it out now
phi1_vec    = mappedX_SNE(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = mappedX_SNE(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = mappedX_SNE(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = mappedX_SNE(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot it!
figure(501)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')


subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')


subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')


subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')

%% Play with options : ManifoldChart
techniques = {'PCA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 'Sammon', 'Isomap', ...
        'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', ...
        'FastMVU', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', ...
        'LLTSA', 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA'};

% Use Autoencoder from DR Toolbox
[mappedX_MC, mapping_MC] = compute_mapping(X, 'ManifoldChart', 4);

% Try it out now
phi1_vec    = mappedX_MC(:,1);
phi1_matrix = reshape(phi1_vec, size(CSN.z_over_s_mesh));
phi2_vec    = mappedX_MC(:,2);
phi2_matrix = reshape(phi2_vec, size(CSN.z_over_s_mesh));
phi3_vec    = mappedX_MC(:,3);
phi3_matrix = reshape(phi3_vec, size(CSN.z_over_s_mesh));
phi4_vec    = mappedX_MC(:,4);
phi4_matrix = reshape(phi4_vec, size(CSN.z_over_s_mesh));

% And plot it!
figure(501)
subplot(221)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi1_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('1^{st} POD Mode')


subplot(222)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi2_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('2^{st} POD Mode')


subplot(223)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi3_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('3^{rd} POD Mode')


subplot(224)
surf(CSN.z_over_s_mesh, CSN.y_over_d_mesh, phi4_matrix)
view(2); shading flat;
axis([min(min(CSN.z_over_s_mesh)) max(max(CSN.z_over_s_mesh)) 0 2]);
xlabel('z/S'); ylabel('y/\delta'); %colorbar;
title('4^{th} POD Mode')
