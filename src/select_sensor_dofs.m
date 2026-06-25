%% Sensor Placement Check With AutoMAC
% Select the measurement DOFs used in the experimental campaign and compute
% AutoMAC for the FE mode shapes at those DOFs.

clear; clc; close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
feaDir = fullfile(rootDir, 'data', 'fea');
figureDir = fullfile(rootDir, 'results', 'figures');

load(fullfile(feaDir, 'modeshapes_full.mat'), 'M');

% Presentation setup: six FE modes and five accelerometer DOFs.
modeColumns = 7:12;
sensorDofs = [12334, 12550, 36526, 10415, 82347];
sensorLabels = {'X1', 'X2', 'X3', 'Y4', 'Z5'};

selectedModeShapes = zeros(numel(sensorDofs), numel(modeColumns));
for i = 1:numel(sensorDofs)
    selectedModeShapes(i, :) = M(sensorDofs(i), modeColumns);
end

AutoMAC = modal_mac(selectedModeShapes, selectedModeShapes);

% Preserve the original variable/file naming used by the old scripts.
modeShape = selectedModeShapes; %#ok<NASGU>
save(fullfile(feaDir, 'modehapes.mat'), 'modeShape');
save(fullfile(feaDir, 'selected_sensor_modes.mat'), ...
    'selectedModeShapes', 'AutoMAC', 'sensorDofs', 'sensorLabels');

figure('Name', 'AutoMAC', 'Color', 'w');
barMAC = bar3(AutoMAC);
for k = 1:length(barMAC)
    zdata = barMAC(k).ZData;
    barMAC(k).CData = zdata;
    barMAC(k).FaceColor = 'interp';
end
colormap(jet);
colorbar;
xlabel('Mode');
ylabel('Mode');
zlabel('AutoMAC');
title('Auto Modal Assurance Criterion');
saveas(gcf, fullfile(figureDir, 'AutoMAC_clean.png'));

fprintf('Selected sensor DOFs:\n');
disp(table(sensorLabels.', sensorDofs.', 'VariableNames', {'Label', 'AbaqusDOF'}));
disp('AutoMAC:');
disp(AutoMAC);

function MAC = modal_mac(phiA, phiB)
    nA = size(phiA, 2);
    nB = size(phiB, 2);
    MAC = zeros(nA, nB);
    for i = 1:nA
        for j = 1:nB
            numerator = abs(phiA(:, i)' * phiB(:, j))^2;
            denominator = (phiA(:, i)' * phiA(:, i)) * (phiB(:, j)' * phiB(:, j));
            MAC(i, j) = numerator / denominator;
        end
    end
end
