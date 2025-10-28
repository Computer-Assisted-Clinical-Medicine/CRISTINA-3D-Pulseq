% -------------------------------------------------------------------------
% CRISTINA Pulseq Implementation
% Author: Valentin Jost (2025)
%
% This code is provided for research and educational use only.
% The author assumes no responsibility for its use or any consequences
% arising from running it on MRI hardware.
%
% Users must verify all safety limits (SAR, gradients, duty cycle, etc.)
% before execution on any scanner.
%
% Provided "as is" without warranty of any kind. Use at your own risk.
%
% Based on: Hösl et al., MRM 2023, DOI: 10.1002/mrm.28284
% -------------------------------------------------------------------------

clear; clc; close all;

gamma = 11.262e6; %MHz/T %Sodium

sys = mr.opts('B0', 3,'MaxGrad', 45*0.8, 'GradUnit', 'mT/m', ...
    'MaxSlew', 0.8*200, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma',gamma);

seq=mr.Sequence(sys);           % Create a new sequence object

%Parameters
rf_dur = 500e-6;             
alpha=90; % flip angle

fov = [240e-3 240e-3 240e-3];     % Define FOV
Nx = 24; Ny = 14; Nz = 8;        % Define FOV and resolution
Tread = 3e-3;
Tpre = 0.8e-3;
Tre = 0.6e-3;

TR = 130e-3; %atm do +10ms to what you want to achieve

% % Define other gradients and ADC events
deltak=1./fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);
gx_re = mr.makeTrapezoid('x',sys,'Area',-gx.area,'Duration',Tre);

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

%Cristina params
Phi0_1 = 90;%degrees
Phi0_2 = 0;%degrees
dPhi = 60; %degrees
nPhi = 360/dPhi;
xi90 = 90; %degrees
xi0 = 0; %degrees
tevo = 10e-3;
tmix = 100e-6;
nTE = 10;

%prepare delays
dTE = gx.flatTime + 2*gx.riseTime + gx_re.flatTime +  2*gx_re.riseTime;
TRrest = (TR - 3*rf_dur - tevo - tmix - nTE*dTE);
delayTR = mr.makeDelay(TRrest);

tevo_delay = mr.makeDelay(tevo);
tmix_delay = mr.makeDelay(tmix);


% Loop over phase encodes and define sequence blocks
% Comment one loop if both phase cycles should be in different file

%Xi0
for iZ=1:Nz
    gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);

    for iY=1:Ny
        gyPre = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);

        for iPhi = 0:nPhi-1
        
            %calc pulses
            Phi = Phi0_1 + dPhi*iPhi;
            [rf1, rf1Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(Phi),'use','excitation','system',sys);
            [rf2, rf2Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(Phi+xi0),'use','other','system',sys);
            [rf3, rf3Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(0),'use','other','system',sys);
            
            % Excitation
            seq.addBlock(rf1);
            seq.addBlock(tevo_delay);
            seq.addBlock(rf2);
            seq.addBlock(tmix_delay);
            seq.addBlock(rf3)

            % Encoding
            seq.addBlock(gxPre,gyPre,gzPre);
            for iTE = 1:nTE
                seq.addBlock(gx,adc);
                seq.addBlock(gx_re);
            end
            seq.addBlock(delayTR)
        end
    end
end

%Xi90
for iZ=1:Nz
    gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);

    for iY=1:Ny
        gyPre = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);

        for iPhi = 0:nPhi-1

            %calc pulses
            Phi = Phi0_1 + dPhi*iPhi;
            [rf1, rf1Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(Phi),'use','excitation','system',sys);
            [rf2, rf2Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(Phi+xi90),'use','other','system',sys);
            [rf3, rf3Delay] = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'PhaseOffset',deg2rad(0),'use','other','system',sys);

            % Excitation
            seq.addBlock(rf1);
            seq.addBlock(tevo_delay);
            seq.addBlock(rf2);
            seq.addBlock(tmix_delay);
            seq.addBlock(rf3)

            % Encoding
            seq.addBlock(gxPre,gyPre,gzPre);
            for iTE = 1:nTE
                seq.addBlock(gx,adc);
                seq.addBlock(gx_re);
            end
            seq.addBlock(delayTR)
        end
    end
end

%% check whether the timing of the sequence is correct
fprintf('Timing check...')
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf(' passed successfully\n');
else
    fprintf(' failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% seq.plot();
return

%% Report
rep = seq.testReport;
fprintf([rep{:}]);
return

%% write sequence
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'CRISTINA3d');
seq.setDefinition('mat', [Nx Ny Nz nPhi nTE]);

% savestr = sprintf('//Cristina3d_xi0.seq');
% savestr = sprintf('//Cristina3d_xi90.seq');
savestr = sprintf('//Cristina3d.seq');
seq.write(savestr);
