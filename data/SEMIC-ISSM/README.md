# Comparison of ISSM-SEMIC Versions

On 26 May 2026, Mathieu implemented SEMIC directly in C++ within ISSM. As a result, ISSM no longer supports the external Fortran version of SEMIC. See the following commits:

* ISSM: Commit 55632d1be  
  https://github.com/ISSMteam/ISSM/commit/55632d1be0870e831f490322ad0ec9476e3c9e9b

* ISSM: Commit 893a94277b2  
  https://github.com/ISSMteam/ISSM/commit/893a94277b269f8f5f641281063100a98b8ad8a9

Here, we compare these two versions of ISSM-SEMIC. First, compile ISSM-SEMIC using:

1. The external Fortran version of SEMIC.
2. The integrated C++ version of SEMIC.

* Install the Fortran and C++ versions of SEMIC, respectively.

```bash
# Build ISSM-893a942 (Fortran SEMIC)
bash install-issm-semic-fortran.sh
```

```bash
# Build ISSM-4da5301 (C++ SEMIC)
bash install-issm-semic-cpp.sh
```

* Make results:

```bash
python GenerGenerateSemicISSM.pyate
```

* Draw figures for comparision using `plot_SEMIC_ISSM_diff.ipynb`.