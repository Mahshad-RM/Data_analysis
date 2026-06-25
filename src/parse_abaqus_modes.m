%% Parse Abaqus Modal Output
% Extract mode-shape vectors from the Abaqus text output used in the
% finite-element modal analysis.
%
% The repository already contains the parsed file:
%   data/fea/modeshapes_full.mat
%
% Run this script only if you want to regenerate it from the Abaqus text
% output.

clear; clc; close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
feaDir = fullfile(rootDir, 'data', 'fea');

candidateFiles = dir(fullfile(feaDir, 'Job-5mm-INP*.txt'));
if isempty(candidateFiles)
    error('No Abaqus text file matching Job-5mm-INP*.txt was found in %s', feaDir);
end

abaqusFile = fullfile(candidateFiles(1).folder, candidateFiles(1).name);
nEigenvalues = 12;

fprintf('Reading Abaqus output:\n  %s\n', abaqusFile);

fileID = fopen(abaqusFile, 'r');
if fileID < 0
    error('Could not open %s', abaqusFile);
end
fileContent = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
lines = fileContent{1};

M = [];

for eigenvalue = 1:nEigenvalues
    startIndex = [];
    endIndex = numel(lines);
    marker = ['E I G E N V A L U E    N U M B E R     ' num2str(eigenvalue)];

    for i = 1:numel(lines)
        if startsWith(lines{i}, marker)
            startIndex = i + 15;
            break;
        end
    end

    if isempty(startIndex)
        warning('Eigenvalue block %d was not found.', eigenvalue);
        continue;
    end

    for j = startIndex + 1:numel(lines)
        if startsWith(lines{j}, 'MAXIMUM')
            endIndex = j - 2;
            break;
        end
    end

    blockLines = lines(startIndex:endIndex);
    blockData = [];

    for k = 1:numel(blockLines)
        values = str2double(strsplit(strtrim(blockLines{k})));
        values = values(~isnan(values));
        if numel(values) >= 4
            blockData(end + 1, :) = values(2:4); %#ok<SAGROW>
        end
    end

    M(:, end + 1) = reshape(blockData.', [], 1); %#ok<SAGROW>
end

outputFile = fullfile(feaDir, 'modeshapes_full.mat');
save(outputFile, 'M');

fprintf('Saved parsed mode shapes to:\n  %s\n', outputFile);
