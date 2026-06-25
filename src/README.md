# MATLAB Scripts

Recommended order:

1. `select_sensor_dofs.m`
   - Loads FE mode shapes.
   - Selects the five experimental sensor DOFs.
   - Computes AutoMAC.

2. `run_modal_analysis.m`
   - Loads processed FRFs from `data/processed_frf`.
   - Converts FRFs to IRFs.
   - Applies Ibrahim Time Domain modal extraction.
   - Compares experimental modes with FE modes using MAC and CoMAC.

Optional scripts:

- `parse_abaqus_modes.m`: regenerates `modeshapes_full.mat` from the Abaqus text export.
- `estimate_frf_from_raw.m`: regenerates processed FRFs if the missing raw impact-test files are later added under `data/raw`.

The old project scripts were interactive and path-dependent. These cleaned scripts use paths relative to the repository root.
