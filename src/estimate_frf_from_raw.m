%% Estimate FRFs From Raw Impact-Hammer Tests
% Optional preprocessing script.
%
% The original raw impact-test files are not present in this repository.
% The processed FRFs used by the modal-analysis script are already saved as:
%   data/processed_frf/P2.mat
%   data/processed_frf/P3.mat
%   data/processed_frf/P4.mat
%
% If raw data are later added, place each impact point in a subfolder:
%   data/raw/P2/*.mat
%   data/raw/P3/*.mat
%   data/raw/P4/*.mat
%
% Each raw .mat file is expected to contain variable Dati, with the force
% channel in column 1 and accelerometer channels in columns 2:end.

clear; clc; close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
rawDir = fullfile(rootDir, 'data', 'raw');
outDir = fullfile(rootDir, 'data', 'processed_frf');

if ~isfolder(rawDir)
    error(['Raw data folder was not found: %s\n' ...
           'The repository contains processed FRFs, so run run_modal_analysis.m instead.'], rawDir);
end

config.fsamp = 10240;
config.hammerSensitivity = 2.361e-3;       % V/N
config.accelerometerSensitivity = 10.2e-3; % V/(m/s^2)
config.triggeredDuration = 10;             % s
config.pretriggerDuration = 0.01;          % s
config.forceWindowAlpha = 0.5;
config.expWindowEndRatio = 0.1;

impactPoints = {'P2', 'P3', 'P4'};
for idx = 1:numel(impactPoints)
    impactPoint = impactPoints{idx};
    inputDir = fullfile(rawDir, impactPoint);
    if ~isfolder(inputDir)
        warning('Skipping %s: folder not found.', inputDir);
        continue;
    end

    [H1, H2, coherence, freq] = estimate_point_frf(inputDir, config); %#ok<ASGLU>
    save(fullfile(outDir, [impactPoint '.mat']), 'H1', 'H2', 'coherence', 'freq');
    fprintf('Saved processed FRFs for %s.\n', impactPoint);
end

function [H1, H2, coherence, freq] = estimate_point_frf(inputDir, config)
    files = dir(fullfile(inputDir, '*.mat'));
    if isempty(files)
        error('No raw .mat files found in %s', inputDir);
    end

    sample = load(fullfile(files(1).folder, files(1).name));
    Dati = sample.Dati;
    nChannels = size(Dati, 2);
    nOutputs = nChannels - 1;
    nTriggered = config.triggeredDuration * config.fsamp;
    nPretrigger = ceil(config.pretriggerDuration * config.fsamp);
    df = 1 / config.triggeredDuration;
    freq = 0:df:config.fsamp / 2;

    Sxx = zeros(nTriggered, nOutputs);
    Sxy = zeros(nTriggered, nOutputs);
    Syx = zeros(nTriggered, nOutputs);
    Syy = zeros(nTriggered, nOutputs);

    for k = 1:numel(files)
        raw = load(fullfile(files(k).folder, files(k).name));
        Dati = raw.Dati;

        force = Dati(:, 1) / config.hammerSensitivity;
        response = Dati(:, 2:end) / config.accelerometerSensitivity;

        triggerLevel = 0.1 * max(abs(force));
        triggerIndex = find(abs(force) >= triggerLevel, 1, 'first') - nPretrigger;
        triggerIndex = max(triggerIndex, 1);

        stopIndex = triggerIndex + nTriggered - 1;
        if stopIndex > numel(force)
            warning('Skipping %s: not enough samples after trigger.', files(k).name);
            continue;
        end

        force = force(triggerIndex:stopIndex);
        response = response(triggerIndex:stopIndex, :);
        [~, peakIndex] = max(abs(force));

        forceWindow = zeros(nTriggered, 1);
        forceWindow(1:min(2 * peakIndex, nTriggered)) = ...
            tukeywin(min(2 * peakIndex, nTriggered), config.forceWindowAlpha);
        force = force .* forceWindow;

        responseWindow = ones(nTriggered, 1);
        nExp = nTriggered - peakIndex + 1;
        responseWindow(peakIndex:end) = exp_window(nExp, config.expWindowEndRatio);
        response = response .* responseWindow;

        inputSpectrum = fft(force) ./ nTriggered;
        outputSpectrum = fft(response) ./ nTriggered;

        Sxx = Sxx + inputSpectrum .* conj(inputSpectrum);
        Syy = Syy + outputSpectrum .* conj(outputSpectrum);
        Sxy = Sxy + outputSpectrum .* conj(inputSpectrum);
        Syx = Syx + inputSpectrum .* conj(outputSpectrum);
    end

    nFiles = numel(files);
    Gxx = one_sided(Sxx, nFiles);
    Gyy = one_sided(Syy, nFiles);
    Gxy = one_sided(Sxy, nFiles);
    Gyx = one_sided(Syx, nFiles);

    H1 = Gxy ./ Gxx;
    H2 = Gyy ./ Gyx;
    coherence = abs(Gxy).^2 ./ (Gxx .* Gyy);
end

function G = one_sided(S, nAverage)
    G = S(1:end / 2 + 1, :) ./ nAverage;
    G(2:end - 1, :) = 2 * G(2:end - 1, :);
end

function win = exp_window(nSamples, endRatio)
    if endRatio == 0
        endRatio = eps;
    end
    tau = -(nSamples / 4) / log(endRatio);
    win = exp(-(1:nSamples).' / tau);
end
