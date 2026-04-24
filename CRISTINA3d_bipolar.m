% -------------------------------------------------------------------------
% CRISTINA Pulseq Implementation with bipolar readout
% Author: Valentin Jost (2025)
%
% Generates a Pulseq .seq file for the CRISTINA triple-quantum sodium MRI
% sequence with a bipolar (alternating-polarity) frequency-encoding readout.
% Odd echoes use the forward gradient (gx), even echoes use the reverse
% gradient (gx_r); all echoes are collected under a single always-on ADC.
% Supports two phase cycles (xi90 and xi0) — either or both can be written
% into the same output file.
%
% This code is provided for research and educational use only.
% Users must verify all safety limits (SAR, gradients, duty cycle, etc.)
% before execution on any scanner. Provided "as is" without warranty.
%
% Based on:
%   Hoesl et al., 2022, https://doi.org/10.1016/j.zemedi.2021.09.001
% -------------------------------------------------------------------------

clear; clc; close all;

%% System Configuration

gamma = 11.262e6;   % Gyromagnetic ratio for 23Na [Hz/T]

% Adapt MaxGrad, MaxSlew, and RF/ADC dead times to your scanner!
sys = mr.opts('B0', 3, 'MaxGrad', 45, 'GradUnit', 'mT/m', ...
              'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 10e-6, 'rfDeadTime', 100e-6, ...
              'adcDeadTime', 10e-6, 'gamma', gamma);

seq = mr.Sequence(sys);

%% Sequence Parameters

rf_dur = 500e-6;                        % RF pulse duration [s] — measure coil ringdown!
alpha  = 90;                            % Flip angle [degrees]
fov    = [240e-3, 240e-3, 240e-3];      % Field of view [m]
Nx = 24;  Ny = 24;  Nz = 24;           % Matrix size [freq, phase, partition]

Tread  = round(3.3e-3 / Nx, 5) * Nx;   % Readout duration [s], integer multiple of dwell time
tdwell = Tread / Nx;                    % ADC dwell time per sample [s]
Tpre   = 0.6e-3;                        % Pre-phaser duration [s]
TR     = 120e-3;                        % Repetition time [s]

% Set xi90=1 and/or xi0=1 to include the respective phase cycle in the output
xi90 = 1;
xi0  = 0;

%% Gradient and ADC Events

deltak = 1 ./ fov;   % k-space step size [1/m]

% Forward and reverse readout trapezoids — RiseTime fixed to one dwell time
% to ensure a uniform flat-top duration for the always-on ADC
gx   = mr.makeTrapezoid('x', sys, 'FlatArea',  Nx*deltak(1), 'FlatTime', Tread, 'RiseTime', tdwell);
gx_r = mr.makeTrapezoid('x', sys, 'FlatArea', -Nx*deltak(1), 'FlatTime', Tread, 'RiseTime', tdwell);

g_dur = mr.calcDuration(gx);   % Duration of one echo lobe [s]

% Pre-phaser: shifts k-space to -kx_max before the echo train
gxPre = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', Tpre);

areaY = ((0:Ny-1) - Ny/2) * deltak(2);   % Phase-encode areas [1/m]
areaZ = ((0:Nz-1) - Nz/2) * deltak(3);   % Partition-encode areas [1/m]

%% CRISTINA Triple-Quantum Filter Parameters

Phi0_xi90 = 90;         % Initial RF phase for xi90 cycle [degrees]
Phi0_xi0  = 0;          % Initial RF phase for xi0  cycle [degrees]
dPhi = 60;              % RF phase increment per cycle step [degrees]
nPhi = 360 / dPhi;      % Number of phase cycle steps

tevo = 10e-3;           % Evolution time between RF1 and RF2 [s]
tmix = 100e-6;          % Mixing time between RF2 and RF3 [s]
nTE  = 25;              % Echo train length

%% Always-On ADC Setup
% A single long ADC window covers the entire echo train. Gradient polarity
% alternates per echo; the ADC remains open throughout. Sample count N is
% rounded up to the nearest adcSamplesDivisor boundary.

gradient_delay = 0;   % Optional offset to account for gradient system delays [s]
adc_length = nTE * g_dur + gradient_delay;

N = floor(adc_length / tdwell);
N = N + (sys.adcSamplesDivisor - mod(N, sys.adcSamplesDivisor));   % Align to divisor

adc = mr.makeAdc(N, 'Dwell', tdwell, 'Delay', tdwell/2, 'system', sys);

%% Delay Preparation

dTE    = mr.calcDuration(gx);

TRrest = TR - 3*(rf_dur + sys.rfRingdownTime + sys.rfDeadTime) ...
            - tevo - tmix - Tpre - nTE*dTE;

delayTR    = mr.makeDelay(TRrest);
tevo_delay = mr.makeDelay(tevo);
tmix_delay = mr.makeDelay(tmix);

% Zero-area placeholder gradient spanning the full echo train window,
% used as the base object for mr.addGradients
ro_spoke_actual = nTE * g_dur + Tpre;
g = mr.makeTrapezoid('x', 'system', sys, 'Duration', ro_spoke_actual, 'Area', 0);

%% Sequence Loops
% Comment out one block below if the two cycles should go into separate files.

% --- Xi90: RF2 phase-shifted by +90° relative to RF1 (TQ coherence pathway) ---
if xi90 == 1
    for iZ = 1:Nz
        gzPre = mr.makeTrapezoid('z', 'Area', areaZ(iZ), 'Duration', Tpre, 'system', sys);

        for iY = 1:Ny
            gyPre = mr.makeTrapezoid('y', 'Area', areaY(iY), 'Duration', Tpre, 'system', sys);

            for iPhi = 0:nPhi-1
                Phi = Phi0_xi90 + dPhi*iPhi;

                [rf1, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(Phi),      'system', sys, 'use', 'excitation');
                [rf2, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(Phi + 90), 'system', sys, 'use', 'other');
                [rf3, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(0),        'system', sys, 'use', 'other');

                % TQ filter: RF1 — tevo — RF2 — tmix — RF3
                seq.addBlock(rf1);
                seq.addBlock(tevo_delay);
                seq.addBlock(rf2);
                seq.addBlock(tmix_delay);
                seq.addBlock(rf3);

                grad_forward = gx;
                grad_reverse = gx_r;

                % Spatial encoding pre-phasers (x/y/z applied simultaneously)
                seq.addBlock(gxPre, gyPre, gzPre);

                % Build composite gradient waveform for the bipolar echo train:
                % odd echoes use grad_forward, even echoes use grad_reverse,
                % each delayed by (echo_idx-1)*g_dur to place them in sequence.
                Grad_array = {g};
                for echo_idx = 1:nTE
                    if mod(echo_idx, 2) == 1
                        grad_forward.delay = (echo_idx-1) * g_dur;
                        Grad_array = [Grad_array, grad_forward];
                    else
                        grad_reverse.delay = (echo_idx-1) * g_dur;
                        Grad_array = [Grad_array, grad_reverse];
                    end
                end

                g_comb = mr.addGradients(Grad_array, sys);
                seq.addBlock(g_comb, adc);
                seq.addBlock(delayTR);
            end
        end
    end
end

% --- Xi0: RF2 has the same phase as RF1 (reference pathway) ---
if xi0 == 1
    for iZ = 1:Nz
        gzPre = mr.makeTrapezoid('z', 'Area', areaZ(iZ), 'Duration', Tpre, 'system', sys);

        for iY = 1:Ny
            gyPre = mr.makeTrapezoid('y', 'Area', areaY(iY), 'Duration', Tpre, 'system', sys);

            for iPhi = 0:nPhi-1
                Phi = Phi0_xi0 + dPhi*iPhi;

                [rf1, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(Phi),     'system', sys, 'use', 'excitation');
                [rf2, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(Phi + 0), 'system', sys, 'use', 'other');
                [rf3, ~] = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, ...
                    'PhaseOffset', deg2rad(0),       'system', sys, 'use', 'other');

                % TQ filter: RF1 — tevo — RF2 — tmix — RF3
                seq.addBlock(rf1);
                seq.addBlock(tevo_delay);
                seq.addBlock(rf2);
                seq.addBlock(tmix_delay);
                seq.addBlock(rf3);

                grad_forward = gx;
                grad_reverse = gx_r;

                % Spatial encoding pre-phasers (x/y/z applied simultaneously)
                seq.addBlock(gxPre, gyPre, gzPre);

                % Build composite gradient waveform for the bipolar echo train
                Grad_array = {g};
                for echo_idx = 1:nTE
                    if mod(echo_idx, 2) == 1
                        grad_forward.delay = (echo_idx-1) * g_dur;
                        Grad_array = [Grad_array, grad_forward];
                    else
                        grad_reverse.delay = (echo_idx-1) * g_dur;
                        Grad_array = [Grad_array, grad_reverse];
                    end
                end

                g_comb = mr.addGradients(Grad_array, sys);
                seq.addBlock(g_comb, adc);
                seq.addBlock(delayTR);
            end
        end
    end
end

%% Timing Check

fprintf('Timing check...');
[ok, error_report] = seq.checkTiming;
if ok
    fprintf(' passed successfully\n');
else
    fprintf(' failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

return

%% Report

rep = seq.testReport;
fprintf([rep{:}]);

return

%% Write Sequence File

seq.setDefinition('FOV',  fov);
seq.setDefinition('Name', 'CRISTINA3d');
seq.setDefinition('mat',  [Nx, Ny, Nz, nPhi, nTE]);
seq.setDefinition('FreqEncDir', grad_forward.channel);

if     xi0 == 1 && xi90 == 0
    savestr = '\\Cristina3d_xi0_bipolar_iso.seq';
elseif xi0 == 0 && xi90 == 1
    savestr = '\\Cristina3d_xi90_bipolar_iso.seq';
elseif xi0 == 1 && xi90 == 1
    savestr = '\\Cristina3d_bipolar_iso.seq';
end

seq.write(savestr);
% safety_check(savestr, PATH_TO_YOUR_.ASC_FILE);