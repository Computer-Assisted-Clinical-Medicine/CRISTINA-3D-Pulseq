% -------------------------------------------------------------------------
% CRISTINA 3D Reconstruction — Bipolar Readout
% Author: Valentin Jost (2025)
%
% Reconstructs single-quantum (SQ) and triple-quantum (TQ) sodium MRI images
% from CRISTINA sequence data acquired with a bipolar readout gradient.
%
% Workflow:
%   1. Load raw k-space data for both phase cycles (xi90, xi0) via mapVBVD
%   2. Read acquisition geometry from the corresponding .seq file
%   3. Deconstruct the always-on ADC stream into a 6-D k-space array
%      [kx, ky, kz, echo, phase_cycle, average]
%   4. Reverse even-echo readout lines to correct bipolar polarity
%   5. Average over repeated acquisitions
%   6. Reconstruct images via 3D iFFT and compute SQ/TQ contrasts via
%      quadrature combination of the two phase cycles
%
% Inputs (adapt file paths before use):
%   - Cristina3d_xi90_bipolar_iso.dat  : raw data, xi90 phase cycle
%   - Cristina3d_xi0_bipolar_iso.dat   : raw data, xi0 phase cycle
%   - Cristina3d_xi0_bipolar_iso.seq   : Pulseq file (geometry/timing source)
%
% Outputs:
%   - SQ : single-quantum image magnitude (echo 5, first time point)
%   - TQ : triple-quantum image magnitude (echo 1, max over echo train)
%   - Figure with central-slice SQ and TQ images
%
% This code is provided for research and educational use only.
%
% Based on:
%   Hoesl et al., 2022, https://doi.org/10.1016/j.zemedi.2021.09.001
% -------------------------------------------------------------------------

clear; clc; close all;

%% Load Raw K-Space Data

% mapVBVD may return a cell array for multi-RAID files; always use the last entry
twix_obj = mapVBVD('\\Cristina3d_xi90_bipolar_iso.dat');
if iscell(twix_obj)
    ksp_Xi90 = twix_obj{end}.image.unsorted();
else
    ksp_Xi90 = twix_obj.image.unsorted();
end

twix_obj = mapVBVD('\\Cristina3d_xi0_bipolar_iso.dat');
if iscell(twix_obj)
    ksp_Xi0 = twix_obj{end}.image.unsorted();
else
    ksp_Xi0 = twix_obj.image.unsorted();
end

clear twix_obj;

%% Read Sequence Parameters

% Load the .seq file to extract FOV, matrix size, and k-space trajectory.
% The xi0 file is used as the geometry reference (both cycles share parameters).
seqpath = '\\Cristina3d_xi0_bipolar_iso.seq';
seq = mr.Sequence();
seq.read(seqpath, 'detectRFuse');

fov               = seq.getDefinition('FOV');
ro_os             = seq.getDefinition('ro');           % Readout oversampling factor
data_size         = seq.getDefinition('mat');          % As stored: [Nx, Ny, Nz, nTE, nPhi]
readout_direction = seq.getDefinition('FreqEncDir');

% Determine number of averages from the ratio of acquired ADC shots to
% expected shots (Ny * Nz * nPhi per average)
avg = size(ksp_Xi0, 3) / (data_size(1) * data_size(3) * data_size(4));

% Reorder and expand to 6-D: [Nx, Ny, Nz, nPhi, nTE, avg]
data_size = [data_size(1), data_size(2), data_size(3), data_size(5), data_size(4), avg];

% Calculate ADC sample times and excitation times from the sequence trajectory
[ktraj_adc, t_adc, ~, ~, t_excitation] = seq.calculateKspacePP();

% Extract echo times relative to the third excitation pulse (RF3), offset
% by 250 µs since a block pulse is used and Pulseq counts from pulse centre.
% Indices 13:24:(24*nTE) pick the centre sample of each echo lobe.
t_echos = t_adc(13:24:(24*data_size(4))) - t_excitation(3) - 250e-6;

%% Reshape ADC Stream into 6-D K-Space Array
% The always-on ADC produces one long sample vector per shot. This loop
% demultiplexes it into [kx, ky, kz, nTE, nPhi, avg], reversing even echoes
% to correct for the alternating readout gradient polarity.

ksp_reshp_Xi90 = zeros(data_size);
ksp_reshp_Xi0  = zeros(data_size);
t              = zeros(data_size(1:5));   % Sample timestamps [kx, ky, kz, nTE, nPhi]

count_adc = 1;   % ADC shot index (increments once per phase-encode + cycle step)

for avg_i = 1:data_size(6)
    for iZ = 1:data_size(3)
        for iY = 1:data_size(2)
            for iPhi = 1:data_size(5)

                % Extract the full ADC readout vector for this shot
                adc_xi90 = ksp_Xi90(:, count_adc);
                adc_xi0  = ksp_Xi0(:,  count_adc);

                count_sample = 1;   % Sample pointer within the current ADC vector

                for iTE = 1:data_size(4)

                    % Slice out Nx samples for this echo
                    samples_xi90 = adc_xi90(count_sample : count_sample + data_size(2) - 1);
                    samples_xi0  = adc_xi0( count_sample : count_sample + data_size(2) - 1);

                    if mod(iTE, 2) == 1
                        % Odd echo: forward gradient — store as acquired
                        ksp_reshp_Xi90(:, iY, iZ, iTE, iPhi, avg_i) = samples_xi90;
                        ksp_reshp_Xi0( :, iY, iZ, iTE, iPhi, avg_i) = samples_xi0;
                    else
                        % Even echo: reverse gradient — flip kx to restore
                        % consistent k-space traversal direction
                        ksp_reshp_Xi90(:, iY, iZ, iTE, iPhi, avg_i) = flip(samples_xi90);
                        ksp_reshp_Xi0( :, iY, iZ, iTE, iPhi, avg_i) = flip(samples_xi0);
                    end

                    % Store ADC sample timestamps on the first average only
                    if avg_i == 1
                        t(:, iY, iZ, iTE, iPhi) = t_adc(count_sample : count_sample + data_size(2) - 1);
                    end

                    % Advance sample pointer: Nx samples + 2 gradient ramp samples per echo
                    count_sample = count_sample + data_size(2) + 2;
                end

                count_adc = count_adc + 1;
            end
        end
    end
end

% Average over repeated acquisitions and remove the singleton average dimension
ksp_reshp_Xi90 = squeeze(mean(ksp_reshp_Xi90, 6));
ksp_reshp_Xi0  = squeeze(mean(ksp_reshp_Xi0,  6));

% 180° k-space rotation to correct for sequence-specific phase accumulation
ksp_reshp_Xi90 = rot90(ksp_reshp_Xi90, 2);
ksp_reshp_Xi0  = rot90(ksp_reshp_Xi0,  2);

%%%OPTIONAL: zero-fill k-space and apply low-pass filters before reconstruction%%%

return

%% Image Reconstruction and Phase-Cycle Combination

fprintf('xi0...');
img_xi0  = ifft3c_new(ksp_reshp_Xi0);
spec_xi0 = fftc_new(img_xi0, 5);     % FFT along echo dimension to get spectral domain

fprintf('xi90...');
img_xi90  = ifft3c_new(ksp_reshp_Xi90);
spec_xi90 = fftc_new(img_xi90, 5);

% Quadrature combination of the two phase cycles yields complex TQ/SQ-weighted images.
% xi0 contributes the real part, xi90 the imaginary part.
img_comb = spec_xi0 + 1i*spec_xi90;

fprintf('done.\n');


%% Advanced method from original paper
% [SQ,TQ,ZQ,...
%     SQ_TE,TQ_TE,ZQ_TE,...
%     B0,imask,imask2] = CRISTINA_getMQC_3D(ksp_reshp_Xi0,ksp_reshp_Xi90,t_echos,60,10e-3, 1);

%% Extract SQ and TQ Contrasts
%%%OPTIONAL: denoise img_comb before coherence extraction%%%
% [img_denoised,~,~,~] = denoise_recursive_tensor(img_comb,[5,5,5],indices={1:3 4 5});

% spectra indices 3 & 5 correspond to single-quantum (SQ)
% spectra index 1 corresponds to the triple-quantum (TQ)
SQ_TE = squeeze(img_comb(:, :, :, :, 5));   % SQ echo train [x, y, z, echo]
TQ_TE = squeeze(img_comb(:, :, :, :, 1));   % TQ echo train [x, y, z, echo]

SQ = abs(squeeze(SQ_TE(:, :, :, 1)));        % SQ image: first echo, magnitude
TQ = abs(squeeze(max(TQ_TE, [], 4)));         % TQ image: max over echo train, magnitude

%% Display Central Slice

slice = floor(data_size(3) / 2);

figure;
subplot(1, 2, 1); imagesc(SQ(:, :, slice)); title('SQ'); axis image; colorbar;
subplot(1, 2, 2); imagesc(TQ(:, :, slice)); title('TQ'); axis image; colorbar;

return

%%%OPTIONAL: voxel-wise fit procedure for TQ/SQ ratio computation%%%