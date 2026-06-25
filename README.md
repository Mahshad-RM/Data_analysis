# Modal Analysis of a Steel Specimen

Data Analysis project, Politecnico di Milano, June 2024.

Authors:

- Pouria Khajehpour
- Mohammad Ali Raheb
- Muhammad Abdullah Saeed Abdullah
- Mahshad Rastegarmoghaddam

Supervisors:

- Prof. Stefano Manzoni
- Dott. Francescantonio Luca

## Project Summary

This project performs modal analysis of a steel specimen using both finite-element and experimental data.

The workflow combines:

- finite-element modal analysis in Abaqus;
- sensor placement assessment with AutoMAC;
- impact-hammer testing with five accelerometers;
- FRF estimation using `H1`, `H2`, and coherence;
- conversion of FRFs to impulse response functions;
- Ibrahim Time Domain modal extraction;
- validation against FE mode shapes using MAC and CoMAC.

## Main Results

The FE model predicted six main natural frequencies:

| Mode | FE model [Hz] | Experimental modal analysis [Hz] |
| --- | ---: | ---: |
| 1 | 831 | 882.921 |
| 2 | 907 | 1012.39 |
| 3 | 1308 | 1330.12 |
| 4 | 1408 | Not clearly identified |
| 5 | 1495 | Not clearly identified |
| 6 | 1907 | 1908.44 |

Modes 4 and 5 were not clearly identified in the experimental FRFs, likely because their resonances had lower amplitudes or were weakly observable from the selected measurement setup.

## Repository Structure

```text
data/
  fea/
  processed_frf/
docs/
presentation/
results/
  figures/
src/
```

## How to Run

Open MATLAB from the repository root or run scripts directly from `src`.

Recommended:

```matlab
run('src/select_sensor_dofs.m')
run('src/run_modal_analysis.m')
```

The script `src/estimate_frf_from_raw.m` is optional. It requires raw impact-test files that are not currently available in this folder. The processed FRFs used in the final analysis are already included as `P2.mat`, `P3.mat`, and `P4.mat`.

## Data Notes

- `data/processed_frf/P2.mat`, `P3.mat`, `P4.mat`: processed FRFs for the three impact points.
- `data/fea/modeshapes_full.mat`: FE mode-shape matrix parsed from Abaqus.
- `data/fea/modehapes.mat`: selected five sensor DOFs used for FE/experimental comparison.
- `data/fea/nat_freqs_full.mat`: saved poles/natural-frequency candidates from the original stabilization analysis.
- `data/fea/Job-5mm-INP(1) (1).txt`: Abaqus text output used for FE mode-shape extraction.

## Publication Cleanup

The original project folder contained draft scripts, interactive path selection, manual thresholding, duplicate files, and generated figures. This cleaned version keeps the presentation, final data, result figures, and reproducible MATLAB scripts with project-relative paths.
