
%%%%% Script to run Example 1: Surface mesh from cardiac contours %%%%%
%
%  File is part of SurFCo package: www.github.com/benvillard/surfco. 
%
%  Please see the following papers for an explanation of the various 
%  parameters: 
%  
%  [1] B. Villard, V. Grau, and E. Zacur, Surface mesh reconstruction from 
%  cardiac MRI contours, J. Imaging, vol. 4(1), no. 16, 2018.
% 
%   [2]  B. Villard, V. Carapella, R. Ariga,  V. Grau, and E. Zacur, 
%   Cardiac Mesh Reconstruction from Sparse, Heterogeneous Contours. 
%   In: Valdés Hernández M., González-Castro V. (Eds.) Medical Image 
%   Understanding and Analysis. MIUA 2017. 
%   Communications in Computer and Information Science, 
%   Vol. 723. Springer, Cham


%% Load Example Data 1 (Surface Mesh from Cardiac Contours)

load('ExampleData_1');

%% Run Mesh Reconstruction Algorithm

[BLID, ULID] = getLids( ExampleData_1 ); % Obtain Arbitrary lids to close mesh
M = SurFCo( ExampleData_1, 'bLID', BLID, 'uLID', ULID, 'plot' ); 

%% Plot Newly Generated Mesh
figure; 
plotMESH( M ); hold on;
arrayfun(@(i) plot3(ExampleData_1{i}(:,1),...      % Plot Original Contours
                    ExampleData_1{i}(:,2),...
                    ExampleData_1{i}(:,3),'r','linewidth',2),...
                    1:numel(ExampleData_1)); 
