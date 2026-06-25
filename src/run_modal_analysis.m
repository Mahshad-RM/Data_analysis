%% Modal Analysis With Ibrahim Time Domain Method
% Reproducible version of the final project workflow.
%
% Inputs:
%   data/processed_frf/P2.mat, P3.mat, P4.mat
%   data/fea/modehapes.mat
%
% Outputs:
%   results/modal_parameters.mat
%   results/figures/*.png

clear; clc; close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
frfDir = fullfile(rootDir, 'data', 'processed_frf');
feaDir = fullfile(rootDir, 'data', 'fea');
resultsDir = fullfile(rootDir, 'results');
figureDir = fullfile(resultsDir, 'figures');

if ~isfolder(figureDir)
    mkdir(figureDir);
end

fsamp = 10240;
dt = 1 / fsamp;
triggeredDuration = 10;
nTriggered = triggeredDuration * fsamp;
analysisTime = 0.06797;
svdRank = 60;
targetFrequencies = [882.921, 1012.39, 1330.12, 1908.44];

[H_f, freq] = load_project_frf(frfDir);
[irf, time] = frf_to_irf(H_f, nTriggered, dt);

plot_frf_overview(freq, H_f, figureDir);
plot_irf_examples(time, irf, figureDir);

Hblocks = collect_irf_blocks(irf);
[~, timeIndex] = min(abs(time - analysisTime));

% The stabilization scan follows the same idea as the original script:
% vary the block-Hankel matrix order and identify stable physical poles.
scanRange = max(1, timeIndex - 200):(timeIndex - 12);
[scanFrequencies, scanDamping] = scan_ibrahim_orders(Hblocks, dt, timeIndex, scanRange, svdRank);

meanDamping = estimate_mean_damping(scanFrequencies, scanDamping, targetFrequencies);
plot_stabilization(scanFrequencies, scanDamping, targetFrequencies, meanDamping, scanRange, figureDir);

chosenBlockColumns = timeIndex - 30;
[chosenPoles, chosenModes] = ibrahim_poles(Hblocks, dt, timeIndex, chosenBlockColumns, svdRank);
[experimentalFrequencies, experimentalModes, selectedPoleIndex] = ...
    select_modes(chosenPoles, chosenModes, targetFrequencies);

[FE_modes, MAC, CoMAC] = compare_with_fea(feaDir, experimentalModes);

save(fullfile(resultsDir, 'modal_parameters.mat'), ...
    'experimentalFrequencies', 'experimentalModes', 'selectedPoleIndex', ...
    'targetFrequencies', 'meanDamping', 'MAC', 'CoMAC', 'FE_modes');

plot_mac(MAC, figureDir);

fprintf('\nIdentified modal frequencies [Hz]:\n');
disp(experimentalFrequencies(:));
fprintf('Mean damping ratios from stable pole bands:\n');
disp(meanDamping(:));
fprintf('MAC matrix between experimental and FE modes:\n');
disp(MAC);
fprintf('CoMAC by selected DOF:\n');
disp(CoMAC(:));

function [H_f, freq] = load_project_frf(frfDir)
    p2 = load(fullfile(frfDir, 'P2.mat'), 'H1', 'freq');
    p3 = load(fullfile(frfDir, 'P3.mat'), 'H1', 'freq');
    p4 = load(fullfile(frfDir, 'P4.mat'), 'H1', 'freq');

    % H1 is stored as frequency x output. Transpose to get output x frequency.
    H_f = [p2.H1.'; p3.H1.'; p4.H1.'];
    freq = p2.freq;
end

function [irf, time] = frf_to_irf(H_f, nTriggered, dt)
    twoSidedFrf = [H_f, fliplr(conj(H_f(:, 2:end)))];
    irf = ifft(twoSidedFrf, [], 2) * nTriggered;
    time = (0:size(irf, 2) - 1) * dt;
end

function Hblocks = collect_irf_blocks(irf)
    % 15 IRFs are ordered as three input locations x five output DOFs.
    nTime = size(irf, 2);
    Hblocks = zeros(5, 3, nTime);
    Hblocks(:, 1, :) = irf(1:5, :);
    Hblocks(:, 2, :) = irf(6:10, :);
    Hblocks(:, 3, :) = irf(11:15, :);
end

function [allFreq, allDamping] = scan_ibrahim_orders(Hblocks, dt, maxTimeIndex, scanRange, svdRank)
    allFreq = cell(numel(scanRange), 1);
    allDamping = cell(numel(scanRange), 1);
    for i = 1:numel(scanRange)
        [poles, ~] = ibrahim_poles(Hblocks, dt, maxTimeIndex, scanRange(i), svdRank);
        freq = abs(imag(poles)) / (2 * pi);
        damping = -real(poles) ./ abs(poles);
        valid = isfinite(freq) & isfinite(damping) & freq > 0 & freq < 3000;
        allFreq{i} = freq(valid);
        allDamping{i} = damping(valid);
    end
end

function [poles, modes] = ibrahim_poles(Hblocks, dt, maxTimeIndex, nBlockColumns, svdRank)
    nBlockRows = maxTimeIndex - nBlockColumns + 1;
    H0 = build_block_hankel(Hblocks, nBlockRows, nBlockColumns, 1);
    H1 = build_block_hankel(Hblocks, nBlockRows, nBlockColumns, 2);

    [U, ~, ~] = svd(H0, 'econ');
    rankUsed = min([svdRank, size(U, 2), size(H0, 1), size(H0, 2)]);
    Utrunc = U(:, 1:rankUsed);

    H0r = Utrunc.' * H0;
    H1r = Utrunc.' * H1;
    [eigVectorsReduced, eigValues] = eig(H1r * pinv(H0r));

    poles = log(diag(eigValues)).' / dt;
    modes = Utrunc * eigVectorsReduced;
end

function H = build_block_hankel(Hblocks, nBlockRows, nBlockColumns, startIndex)
    nOutputs = size(Hblocks, 1);
    nInputs = size(Hblocks, 2);
    H = zeros(nOutputs * nBlockRows, nInputs * nBlockColumns);

    for row = 1:nBlockRows
        rowIdx = (row - 1) * nOutputs + (1:nOutputs);
        for col = 1:nBlockColumns
            colIdx = (col - 1) * nInputs + (1:nInputs);
            H(rowIdx, colIdx) = Hblocks(:, :, startIndex + row + col - 2);
        end
    end
end

function meanDamping = estimate_mean_damping(scanFrequencies, scanDamping, targetFrequencies)
    meanDamping = nan(size(targetFrequencies));
    for i = 1:numel(targetFrequencies)
        freqBand = 0.001 * targetFrequencies(i);
        dampingValues = [];
        for k = 1:numel(scanFrequencies)
            idx = abs(scanFrequencies{k} - targetFrequencies(i)) <= freqBand;
            dampingValues = [dampingValues; scanDamping{k}(idx).']; %#ok<AGROW>
        end
        dampingValues = dampingValues(dampingValues > 0 & dampingValues < 0.02);
        meanDamping(i) = mean(dampingValues, 'omitnan');
    end
end

function [frequencies, modeShapes, selectedIndex] = select_modes(poles, modes, targetFrequencies)
    poleFrequencies = abs(imag(poles)) / (2 * pi);
    frequencies = zeros(size(targetFrequencies));
    selectedIndex = zeros(size(targetFrequencies));
    modeShapes = zeros(5, numel(targetFrequencies));

    for i = 1:numel(targetFrequencies)
        [~, idx] = min(abs(poleFrequencies - targetFrequencies(i)));
        selectedIndex(i) = idx;
        frequencies(i) = poleFrequencies(idx);
        modeShapes(:, i) = real(modes(1:5, idx));
    end

    modeShapes = normalize_columns(modeShapes);
end

function [FE_modes, MAC, CoMAC] = compare_with_fea(feaDir, experimentalModes)
    FE = load(fullfile(feaDir, 'modehapes.mat'), 'modeShape');

    % Sign convention and modal columns used in the original final script:
    % experimental modes 1,2,3,6 were identified.
    FE_modes_all = [FE.modeShape(1:3, :); -FE.modeShape(4, :); FE.modeShape(5, :)];
    FE_modes = FE_modes_all(:, [1, 2, 3, 6]);
    FE_modes = normalize_columns(FE_modes);

    MAC = modal_mac(experimentalModes, FE_modes);
    CoMAC = coordinate_mac(experimentalModes, FE_modes);
end

function normalized = normalize_columns(values)
    normalized = values;
    scale = max(abs(values), [], 1);
    scale(scale == 0) = 1;
    for i = 1:numel(scale)
        normalized(:, i) = values(:, i) ./ scale(i);
    end
end

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

function CoMAC = coordinate_mac(phiA, phiB)
    nDof = size(phiA, 1);
    CoMAC = zeros(nDof, 1);
    for dof = 1:nDof
        numerator = sum(phiA(dof, :) .* phiB(dof, :))^2;
        denominator = sum(phiA(dof, :).^2) * sum(phiB(dof, :).^2);
        CoMAC(dof) = numerator / denominator;
    end
end

function plot_frf_overview(freq, H_f, figureDir)
    avgFrf = mean(abs(H_f), 1);
    figure('Name', 'Average FRF', 'Color', 'w');
    semilogy(freq, avgFrf, 'LineWidth', 1.5);
    grid on;
    xlim([0, 3000]);
    xlabel('Frequency [Hz]');
    ylabel('Average |H_1| [m/s^2/N]');
    title('Average Frequency Response Function');
    saveas(gcf, fullfile(figureDir, 'average_frf_clean.png'));
end

function plot_irf_examples(time, irf, figureDir)
    figure('Name', 'Impulse Response Examples', 'Color', 'w');
    exampleIrfs = [1, 8, 15];
    for k = 1:numel(exampleIrfs)
        subplot(numel(exampleIrfs), 1, k);
        irfIndex = exampleIrfs(k);
        plot(time, real(irf(irfIndex, :)), 'LineWidth', 1.0);
        grid on;
        xlim([0, 0.12]);
        xlabel('Time [s]');
        ylabel('IRF');
        title(sprintf('IRF %d', irfIndex));
    end
    saveas(gcf, fullfile(figureDir, 'irf_examples_clean.png'));
end

function plot_stabilization(scanFrequencies, scanDamping, targetFrequencies, meanDamping, scanRange, figureDir)
    figure('Name', 'Stabilization Diagram', 'Color', 'w');
    hold on; grid on;
    for k = 1:numel(scanFrequencies)
        order = scanRange(k);
        freq = scanFrequencies{k};
        damping = scanDamping{k};
        plot(freq, order * ones(size(freq)), '+k', 'LineStyle', 'none');

        for i = 1:numel(targetFrequencies)
            inFreqBand = abs(freq - targetFrequencies(i)) <= 0.001 * targetFrequencies(i);
            inDampingBand = abs(damping - meanDamping(i)) <= 0.05 * meanDamping(i);
            stable = inFreqBand & inDampingBand;
            plot(freq(stable), order * ones(size(freq(stable))), 'ro', 'LineStyle', 'none');
        end
    end
    xlim([0, 3000]);
    xlabel('Frequency [Hz]');
    ylabel('Model order parameter');
    title('Stabilization Diagram');
    legend({'Poles', 'Stable physical poles'}, 'Location', 'best');
    saveas(gcf, fullfile(figureDir, 'stabilization_clean.png'));
end

function plot_mac(MAC, figureDir)
    figure('Name', 'MAC', 'Color', 'w');
    imagesc(MAC);
    axis equal tight;
    colorbar;
    xlabel('FE mode');
    ylabel('Experimental mode');
    title('MAC: Experimental vs FE Modes');
    saveas(gcf, fullfile(figureDir, 'mac_clean.png'));
end
