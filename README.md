# CRISTINA — 3D Sodium Multi-Quantum Coherence MRI (Pulseq)

This repository is part of an abstract submitted to ISMRM 2026. It provides a
ready-to-compile open-source implementation of the 3D sodium multi-quantum
coherence (MQC) sequence **CRISTINA** (*Cartesian imaging of single and
triple-quantum sodium*), implemented in MATLAB using the
[Pulseq](https://pulseq.github.io/) framework, together with a corresponding
reconstruction script.

> **Note on monopolar readout:** The monopolar readout showed large signal
> deviations across the echo dimension. This may be caused by gradient or ADC
> delays that we could not resolve or measure appropriately. We therefore
> adapted the implementation to use a bipolar readout. If you have encountered
> the same issue or have ideas on how to resolve it, please feel free to open
> an issue or reach out — we welcome further discussion.

> **Note on bipolar readout:** The bipolar implementation exhibits a residual
> artefact: a bright line appears at the centre of the image along the
> frequency-encoding direction. This is a known issue in bipolar multi-echo
> acquisitions, likely caused by a DC offset or a k-space centre inconsistency
> between odd (forward) and even (reverse) echoes. It has not yet been fully
> resolved. If you have suggestions or a working correction, please open an
> issue or get in touch.
---

## Overview

CRISTINA is a three-pulse triple-quantum filter sequence that simultaneously
acquires single-quantum (SQ) and triple-quantum (TQ) filtered ²³Na signals via
two complementary phase cycles (xi90 and xi0). This implementation supports:

- **Monopolar readout** — separate rewinder gradient after each echo
- **Bipolar readout** — alternating-polarity gradients with a single always-on ADC window
- **3D Cartesian k-space** sampling with configurable matrix size and echo train length
- Flexible phase-cycle selection: xi90 only, xi0 only, or both in a single file

---

## Dependencies

| Dependency | Purpose |
|---|---|
| [Pulseq for MATLAB](https://github.com/pulseq/pulseq) | Sequence design and `.seq` file export |
| [mapVBVD](https://github.com/CIC-methods/FID-A) | Raw data loading (Siemens TWIX format) |
| [MQC_Imaging_CRISTINA](https://github.com/MHoesl/MQC_Imaging_CRISTINA) | Original CRISTINA reconstruction framework by Hoesl et al. |
| [Tensor-MP-PCA](https://github.com/Neurophysics-CFIN/Tensor-MP-PCA) | Optional k-space denoising (Olesen et al., MRM 2022) |

---

## Usage

### 1. Generate the sequence file

Open `Cristina3d_monopolar.m` or `Cristina3d_bipolar.m` and adapt:

- **System limits** (`MaxGrad`, `MaxSlew`, RF dead times) to your scanner
- **Sequence parameters** (`fov`, `Nx/Ny/Nz`, `TR`, `nTE`, `tevo`, `tmix`)
- **Phase cycle flags** (`xi90`, `xi0`)
- **Output path** (`savestr`)

Run the script. It will perform a Pulseq timing check and, if passed, write
a `.seq` file to the specified path.

> ⚠️ **Safety notice:** Always verify SAR, gradient amplitude, slew rate, and
> duty cycle limits on your specific scanner before executing any sequence on
> hardware. This code is provided "as is" without warranty of any kind.

### 2. Reconstruct

Open `Recon_Cristina3d.m` and adapt the file paths to your `.dat` and `.seq`
files, then run section by section:

1. **Load** raw k-space data via `mapVBVD`
2. **Reshape** the always-on ADC stream; even echoes are flipped to correct
   bipolar polarity
3. **Reconstruct** via 3D iFFT and quadrature combination of the two phase
   cycles
4. **Visualise** the central slice of the SQ and TQ images

## References

1. **Hoesl, M.A.U. et al. (2020).**
   Efficient ²³Na triple-quantum signal imaging on clinical scanners:
   Cartesian imaging of single and triple-quantum ²³Na (CRISTINA).
   *Magnetic Resonance in Medicine*, 84(5), 2412–2428.
   https://doi.org/10.1002/mrm.28284

2. **Hoesl, M.A.U. et al. (2022).**
   Volumetric ²³Na single and triple-quantum imaging at 7T: 3D-CRISTINA.
   *Zeitschrift für Medizinische Physik*, 32(3), 310–323.
   https://doi.org/10.1016/j.zemedi.2021.09.001

3. **Olesen, J.L. et al. (2022).**
   Tensor-MP-PCA denoising for MRI (optional post-processing).
   *Magnetic Resonance in Medicine*, 88(6), 2578–2595.
   https://doi.org/10.1002/mrm.29478
   — Implementation: https://github.com/Neurophysics-CFIN/Tensor-MP-PCA

---

## License

Copyright (C) 2025 Valentin Jost

This program is free software: you can redistribute it and/or modify it under
the terms of the **GNU General Public License version 3** as published by the
Free Software Foundation.

This program is distributed in the hope that it will be useful, but **WITHOUT
ANY WARRANTY**; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html)
for more details.

**Note on dependencies:** [Pulseq](https://github.com/pulseq/pulseq) is
licensed under MIT, which is compatible with GPL-3.0. Any derivative works that
incorporate or modify code from this repository must also be released under
GPL-3.0.
