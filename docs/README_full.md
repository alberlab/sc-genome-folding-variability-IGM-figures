> **Note (GitHub copy, 2026-11-13).** This is the detailed README of the authors' working package
> (`Figures_final_2026-11`). In this repository the layout differs as follows: input tables shared by several figures live
> once in the top-level `data/` (`data/locale/v2_ab4r/`, `data/distCT/e2c/` incl. the per-locus `.npz` and
> `SDrad_vs_SDdistCT_4stage_perbin.tsv`, `data/radial/`, `data/schic/`); `<Figure>/data/` keeps only tables unique to that
> figure; checksums are in `data/checksums.md5`; only the default locale (v2_ab4r) and dist-to-CT normalisation (e2c) are
> shipped; `run.sh` uses `$PYTHON` or the `python` of the active environment (`environment.yml`), not per-figure conda envs
> (their exact records are in `envs/authors_osx-arm64/`); outputs, `data/structures/`, `.mrc` and `.cxs` files are not in
> git (see `fetch_data.sh`). Paths such as `~/plots/...` or `~/Figures_reproducible/...` below refer to the authors' machine.

# Model (IGM) panels of the manuscript figures — reproducible package (2026-11)

This folder rebuilds every **model (IGM) panel** of the current manuscript figures:
`sc-project manuscripts main figures.pdf` and `sc-project manuscripts ext figures.pdf` (version of 2026-11-11). Panels use
the **new panel numbering** of those PDFs. All code and input tables are inside this folder; there are no symlinks and no
HSS files.

## Quick start

```bash
bash run_all.sh                 # every plotted panel (8 figure folders, ~45 s); each folder runs in its own conda env
bash Fig2/run.sh                # one figure only
PYTHON=/path/to/python bash run_all.sh   # one interpreter for everything (needs numpy pandas scipy matplotlib seaborn scikit-learn openpyxl)
bash Fig5/run_chimerax.sh       # ChimeraX renders (not in run_all.sh): Fig5a (SD) + Fig5b (entropy)
bash ExtFig2/run_chimerax.sh    # ExtFig2c locale renders (ES)
bash Fig4/run_chimerax.sh       # Fig4a/e/i: re-render the PNGs from the saved .cxs sessions into Fig4/output/figures/rerender/

# another nuclear-locale version (e.g. NL-act / NL-rep split with a 2.5 r A/B neighbourhood):
bash set_locale_version.sh v2_ab2p5r /Users/yewang/plots/locale_v2_ab2p5r/data/locale_v2_ab2p5r   # copies locale_<cell>.tsv
LOCALE_VERSION=v2_ab2p5r bash run_all.sh                 # -> <Figure>/output_v2_ab2p5r/{figures,source_data}

# another dist-to-CT normalisation (signed distance / bead radius r, or / R_eq of the core hull; e2c = default):
DISTCT_NORM=r   bash run_all.sh                          # -> <Figure>/output_distCT-r/   (tables in data/distCT/r/)
DISTCT_NORM=Req bash run_all.sh                          # -> <Figure>/output_distCT-Req/
LOCALE_VERSION=v2_ab2p5r DISTCT_NORM=Req bash run_all.sh # -> <Figure>/output_v2_ab2p5r_distCT-Req/
```

Each figure folder has the same layout:

```
<Figure>/code/            one plot_<panel>.py per panel (+ the figure's *_common.py; ChimeraX scripts in code/chimerax/)
<Figure>/data/            input tables (real copies) + checksums.md5  (check: cd data && md5sum -c checksums.md5)
<Figure>/data/locale/<version>/<cell>_occ06_locale_N565.tsv   nuclear-locale tables, one folder per locale version
<Figure>/output/figures/  <panel>_<description>.{png,pdf}; ChimeraX renders in sub-folders (Fig4a, Fig5a, ExtFig2c, ...)
<Figure>/output/source_data/  the numbers behind each panel
<Figure>/run.sh           runs the panel scripts of this figure
```

Conda envs (`envs/<env>.yml` = `conda env export --from-history`, `envs/<env>.lock.txt` = full `conda list --export`):
Fig2 `fig2-h1-model`, Fig3 `fig3-4stage-model`, Fig4 `fig4-4stage-model`, Fig5 `fig5-4stage-model`, ExtFig2 `extfig2-model`,
ExtFig4 `extfig4-model`, ExtFig5 `extfig5-model`, ExtFig6 `extfig6-model` (all Python 3.13, pandas 3.0.5,
matplotlib 3.11.1, scipy 1.18.0). `run.sh` uses `$PYTHON`, else `<conda base>/envs/<env>/bin/python`, else `python`.
ChimeraX: UCSF ChimeraX 1.11.1 (auto-detected in /Applications, or set `CHIMERAX=`).

## Conventions (same as the PDFs)

- **Stages / models:** ES, DE, HB, HLC = model ids `h1`, `DE`, `HB`, `mHEP` (IGM occupancy 0.6, 200 kb, 1000 diploid
  structures each). iHLC has no model; panels that show an iHLC column in the PDF leave it blank ("Not count").
- **Bins:** autosomal, non-centromeric 200 kb bins (13,315 per model); 1 Mb values = mean of these 200 kb bins.
- **dist to CT:** signed distance to the chromosome-territory core (convex hull of the top 30 % strict-cis-cisratio beads of
  each chromosome copy, contacts within 2.5 bead radii), normalised by the core centroid-to-surface distance; SD over
  2 copies x 1000 structures per bin (`*_occ06_distCT_top30_2p5r_stats.tsv`). **SD only** — every IQR panel, track, row and
  loop iteration of the earlier code was removed (tables that carry an `iqr` column keep it as data).
- **Locale = ENV v2** (everywhere): speckle 5 r, nucleolus 6 r; env seeds = outermost 20 % of **all** beads (incl.
  centromere / gap beads), env = non-gap seeds + non-gap beads within 2.5 r of a seed; others NL-act / NL-rep by A/B ratio.
  NL-act / NL-rep use the A/B ratio within radius_factor 2 (= 4 r). Source:
  `~/plots/occ06_frank_env_radial_QC/env/data/locale_seed20all_2p5r/locale_<cell>.tsv`, copied as
  `data/locale/v2_ab4r/<cell>_occ06_locale_N565.tsv` (file name kept from the earlier packages; content = ENV v2).
- **Locale version switch:** every `*_common.py` (and `ExtFig2/code/plot_extfig2d.py`, which reads `Fig5/data` directly)
  reads `data/locale/${LOCALE_VERSION:-v2_ab4r}/`. The default writes to `output/`; any other version writes figures and
  source data to `output_<LOCALE_VERSION>/`, so versions never overwrite each other. Install a version with
  `bash set_locale_version.sh <version> <src_dir>`: it checks that `<src_dir>/locale_{h1,DE,HB,mHEP}.tsv` have the v2_ab4r
  header and 15,453 rows, copies them (real files) to Fig2 / ExtFig2 (h1) and Fig3 / Fig4 / Fig5 / ExtFig6 (all four;
  ExtFig2d / 2e and ExtFig6c read `Fig5/data`), and rebuilds every `data/checksums.md5` (no arguments = only rebuild
  checksums). Fig4 uses only the `label` (centromere) column, which is identical between versions; ExtFig4 / ExtFig5 read
  no locale table and always write `output/`. The ChimeraX renders (Fig5b, ExtFig2c) are fixed to v2_ab4r; new renders for
  another version need new cmm/pdb from the cluster (`provenance/`).
- **dist-to-CT normalisation switch (2026-11-12):** `DISTCT_NORM` = `e2c` (default; the shipped
  `data/<cell>_occ06_distCT_top30_2p5r_stats.tsv`, signed distance / core centroid-to-surface distance) | `r` (signed
  distance in nm / the model's bead radius: ES 94.91, DE 106.23, HB 99.03, HLC 108.07 nm) | `Req` (signed distance /
  (3 V_hull / 4π)^(1/3) of the same core hull, per structure and chromosome copy). Same CT core, same 2 x 1000 pooling,
  same table format. Non-default tables live in `data/distCT/<norm>/` (Fig2 / ExtFig2: h1 table + the per-locus npz;
  ExtFig2 also `SDrad_vs_SDdistCT_4stage_perbin.tsv` with only the `sd` column replaced; Fig3 / Fig4 / Fig5 / ExtFig4 /
  ExtFig6: all four tables; ExtFig2d / 2e and ExtFig6c read `Fig5/data`), installed by
  `bash set_distct_norm.sh <r|Req> <src_dir>` (src = `/Users/yewang/plots/distCT_normalisation/data/install/<norm>`; it also
  rebuilds the checksums). Non-default runs write `output_distCT-<norm>/` (after the locale suffix when both are set),
  label axes "… (dist to CT / r)" / "… / R_eq", and multiply the few fixed dist-to-CT axis limits (Fig2c left, Fig2h /
  ExtFig2h histogram range, Fig3d x range, Fig3e y range) by `CT_SCALE` = median over cells of median SD(norm) / median
  SD(e2c) (r 4.79, Req 0.480); fixed histogram y limits become automatic. Not switched: ExtFig5 (ablation / envswap
  models, e2c only) and the Fig5a / 5b ChimeraX renders + colour bars (skipped for r / Req). Computation, gates and the
  comparison of the three normalisations: `/Users/yewang/plots/distCT_normalisation/README.md`.
- **Dominant locale** = argmax of the five exclusive class frequencies (spe, nuc, env, NL-act, NL-rep); Fig5e is the 5-class
  Sankey. Locale entropy = Shannon entropy (bits) of the five frequencies.
- **Colours:** spe #8b1e3f, NL-act #e57373, nuc #f4a261, NL-rep #9b59b6, env #7fa9c9.

## Panel table

Status: **=ref** = source data md5-identical to the reference run **and** PNG pixel-identical; **subset** = same rows as the
reference minus the IQR rows; **new** = new script / layout (checked against the PDF). References: PK =
`~/Figures_reproducible/<pkg>/output`, V2 = `~/plots/Figures_updated_ENV_locale/ENV_v2`, NH4 = 4-stage runs in
`~/plots/Figures_noHB_ENVv2/{Fig3,ExtFig6}/output_withHB` (NH = `~/plots/Figures_noHB_ENVv2/Fig2/output`).

| Panel | Content | Script (`<Figure>/code/`) | Output (`output/figures/`) | Source / reference | Status |
|---|---|---|---|---|---|
| Fig2b | chr11 1 Mb tracks: bulk Eig1, cis/total; plate / droplet Mean & SD(cisratio); model Mean & SD(dist to CT) | `plot_fig2b.py` | `Fig2b_chr11_tracks` | PK Fig2 `plot_fig2b.py` (IQR track removed) | table =ref; figure new (8 tracks) |
| Fig2c | Mean vs Eig1 chr2 (colour SD) / SD vs Eig1 / SD vs Mean (colour Eig1), 1 Mb | `plot_fig2c.py` | `Fig2c_{left,middle,right}_*` | PK Fig2 (SD set only) | =ref, stats subset |
| Fig2d | entropy by dominant locale (box order spe, NL-act, NL-rep, env) / entropy vs SD chr2 / entropy vs Eig1 200 kb | `plot_fig2d.py` | `Fig2d_{left,middle,right}_*` (+ `_alt_`: genome-wide SD, 1 Mb Eig1) | NH Fig2 `plot_fig2d.py` (v2) | =ref, stats subset |
| Fig2e | chr11 1 Mb tracks: locale entropy, RAD, SD(RAD) | `plot_fig2e.py` | `Fig2e_chr11_tracks` | PK Fig2 + v2 locale; ref V2 Fig2_output | =ref |
| Fig2h | model chr11 5 loci inside/outside bars + histogram | `plot_fig2h.py` | `Fig2h_chr11_5loci_bar_hist` | PK Fig2 | =ref |
| Fig3c | IGM rows: Mean & SD(dist to CT), chr2 + chr11, ES DE HB (iHLC blank) HLC | `plot_fig3c.py` | `Fig3c_IGM_all_stages_*` (+ one file per stage) | PK Fig3 (IQR row removed) | table =ref; figure new (2 rows) |
| Fig3d | IGM columns: SD(dist to CT) histogram; entropy histogram + GMM (thick lines lw 8 / 6.5) | `plot_fig3d.py` | `Fig3d_IGM_SDdistCT_hist`, `Fig3d_IGM_entropy_hist_GMM` | NH Fig3 `plot_fig3d.py`, 4 stages; ref NH4 | =ref, summary subset |
| Fig3e | IGM row: Mean(dist to CT) vs own-stage Eig1, chr2, colour SD (1 x 4; the PDF adds the blank iHLC column) | `plot_fig3e.py` | `Fig3e_IGM_Mean_vs_Eig1_chr2_colorbySD` | PK Fig3 | =ref, stats subset |
| Fig4a / 4e / 4i | ChimeraX: nuclear shape half view / chr11 CT / Shannon index | `chimerax/fig4_render_sessions.py` (`run_chimerax.sh`) | `Fig4a/`, `Fig4e/`, `Fig4i/` (PNG + .cxs; Fig4i also the envelope maps `*_nucleus_shell_pores.mrc`) | PK Fig4 renders, copied; cmm/pdb in `data/structures/*.tar.gz` | copied |
| Fig4b | mean SD(dist to CT) per radial decile | `plot_fig4b.py` | `Fig4b_SDdistCT_vs_radial_decile` | PK Fig4 | =ref, shells subset |
| Fig4f | chr11 CT volume ratio / compactness (packing density) violins | `plot_fig4f.py` | `Fig4f_chr11_{CT_volume_ratio,compactness}` | PK Fig4 | =ref |
| Fig4j | mean Shannon index per radial decile | `plot_fig4j.py` | `Fig4j_Shannon_vs_radial_decile` | PK Fig4 | =ref |
| Fig5a | ChimeraX half-sphere, SD(dist to CT) + speckle centroids | `chimerax/fig5ab_*.py` (`run_chimerax.sh SD`); colour bar `plot_fig5ab_colorbars.py` | `Fig5a/` | PK Fig5 SD renders (= V2 FINAL, md5-identical) | copied; colour bar =ref |
| Fig5b | ChimeraX half-sphere, locale entropy (v2) | `chimerax/fig5ab_*.py` (`run_chimerax.sh entropy`) | `Fig5b/` | V2 FINAL `Fig5b_renders` | copied; colour bar =ref |
| Fig5c | composition of all bead instances: spe, NL-act, nuc, NL-rep, env, other | `plot_fig5ce_per_locale.py` | `Fig5c_locale_composition_stacked` | `~/plots/.../new_panels/code/plot_fig5ce_per_locale.py`; ref V2 | =ref |
| Fig5d | entropy vs SD, chr2 / chr11 x 4 stages | `plot_fig5d.py` | `Fig5d_entropy_vs_SDdistCT_chr2_chr11` | `plot_new_panels.py` (Fig5d part); ref V2 | =ref |
| Fig5e | 5-class dominant-locale Sankey ES → DE → HB → HLC | `plot_fig5ce_per_locale.py` | `Fig5e_dominant_locale_sankey` | as Fig5c; ref V2 | =ref |
| ExtFig2a | model SD(dist to CT) vs droplet / plate SD(cisratio), 1 Mb density | `plot_extfig2a.py` | `ExtFig2a_model_vs_schic_SDcisratio` | PK ExtFig2 (2 x 2 → 1 x 2, SD only, droplet first) | new layout; numbers =ref |
| ExtFig2b left | SD(dist to CT) vs RAD, ES density | `plot_extfig2b.py` | `ExtFig2b_left_SDdistCT_vs_RAD` | PK ExtFig2 (SD panel only) | new layout; numbers =ref |
| ExtFig2b right | SD(RAD) vs SD(dist to CT), ES density | `plot_extfig2b_right.py` | `ExtFig2b_right_SDRAD_vs_SDdistCT_ES` | `~/plots/sdrad_vs_sddistct_4stage/code/plot_sdrad_vs_sddistct_density.py` | PNG pixel-identical; r = 0.3727 = ref |
| ExtFig2c | ES locale renders (5 classes) + 3 example loci locale-frequency bars | `plot_extfig2c.py`; `chimerax/extfig2c_*.py` (`run_chimerax.sh`) | `ExtFig2c_example_loci_locale_freq`, `ExtFig2c/` | PK ExtFig2 + v2; renders V2 FINAL (`extfig2c_seed20all_2p5r`) | =ref (V2) |
| ExtFig2d | ES all autosomes, entropy vs SD: active (spe + NL-act dominant) / repressive (env + NL-rep + nuc), opacity = strength | `plot_extfig2d.py` (data `../Fig5/data`) | `ExtFig2d_ES_entropy_vs_SDdistCT_active_repressive` | `~/plots/entropy_CV_vs_SDCT_diagnosis/code/plot_entropy_vs_IQR_by_locale_strength.py sd 2class` (ES column, 1 x 2) | new; v1 numbers reproduced exactly, v2 below |
| ExtFig2e | ES chr11 / chr20 entropy vs SD | `plot_extfig2e.py` (data `../Fig5/data`) | `ExtFig2e_ES_entropy_vs_SDdistCT_chr11_chr20` | `plot_new_panels.py` (old name ExtFig2d); ref V2 | =ref |
| ExtFig2h | model chr2 3 loci inside/outside bars + histogram | `plot_extfig2h.py` | `ExtFig2h_chr2_3loci_bar_hist` | PK ExtFig2 `plot_extfig2e.py` (old name 2e) | =ref |
| ExtFig4d | IGM row: chr11 5 loci, Mean & SD(dist to CT) across stages | `plot_extfig4.py` | `ExtFig4d_IGM_chr11_5loci_Mean_SDdistCT` | PK ExtFig4 (SD figure) | =ref |
| ExtFig5b | model chr2 CT volume ratio / compactness (packing density) violins | `plot_extfig5b.py` | `ExtFig5b_chr2_{CT_volume_ratio,compactness}` | PK ExtFig5 | =ref |
| ExtFig6a | top: SD(dist to CT) vs RAD; bottom: entropy vs RAD; 4 stages, colour = dominant locale | `plot_extfig6ab.py` | `ExtFig6a_{top_SDdistCT,bottom_entropy}_vs_RAD` | PK ExtFig6 + v2; ref V2 | =ref, stats subset |
| ExtFig6b | % of top-10 % SD(dist to CT) bins that are spe-dominant | `plot_extfig6_spe_topvar.py` | `ExtFig6b_top10_SDdistCT_fraction_spe_dominant` (+ `_alt_`: % of spe bins in top 10 %) | PK ExtFig6 + v2; ref V2 | =ref, table subset |
| ExtFig6c | entropy vs SD, all autosomes / chr20, 4 stages | `plot_extfig6c.py` (data `../Fig5/data`) | `ExtFig6c_entropy_vs_SDdistCT_merged_chr20` | `plot_new_panels.py` (ExtFig6c part); ref V2 | =ref |
| ExtFig6d | top: SD(dist to CT) boxes; bottom: entropy boxes + two-sided Mann–Whitney p (each stage vs ES); spe, NL-act, env, NL-rep | `plot_extfig6d.py` | `ExtFig6d_box_SD_entropy_by_dominant_locale_MWU` | new: `draw()` SD row + `draw_entropy_tests()` (4 locales) of NH ExtFig6 `plot_extfig6d.py` on one 2 x 4 grid | new; box table = NH4 (same rows), test table md5 =ref NH4 |
| ExtFig6e | entropy vs SD 4 x 4 (rows spe, env, NL-act, NL-rep), solid colours + legend | `plot_extfig6e_solid.py` | `ExtFig6e_entropy_vs_SD_by_dominant_locale_solid` | NH ExtFig6, 4 stages; ref NH4 | =ref, stats subset |

**Cross-folder reads:** ExtFig2d, ExtFig2e and ExtFig6c read the four-stage tables in `Fig5/data` (ExtFig2e / ExtFig6c through
`Fig5/code/fig5_common.py`, unchanged from the script they were split from). The ChimeraX scenes of Fig5 and ExtFig2c open
the envelope maps in `Fig4/output/figures/Fig4i/`. Keep the folder structure together.

## ChimeraX renders

| Panel | Inputs (`data/structures/`, unpack with `extract_structures.sh`) | Renderer | Shipped PNG |
|---|---|---|---|
| Fig4a / 4e / 4i | `<cell>_occ06_cmm_pdb_{by_chromosome,by_shannon}.tar.gz` (structures 82, 399, 427, 602, 833) | sessions `output/figures/Fig4{a,e,i}/*.cxs`; `run_chimerax.sh` re-saves them as 2000 px PNG into `output/figures/rerender/` (scenes were composed interactively, the session is the recipe) | PK Fig4 |
| Fig5a | `<cell>_occ06_fig5_SD_entropy_ENVv2_cmm_pdb.tar.gz` (`cmm_pdb_SD`, `speckles`) | `code/chimerax/fig5ab_all.py` → `fig5ab_scene.py`, `fig5_speckle_legend.py` | PK Fig5a SD (md5 = V2 FINAL) |
| Fig5b | same tarball (`cmm_pdb_entropy`, B-factor = v2 locale entropy) | same, metric `entropy` | V2 FINAL `Fig5b_renders` |
| ExtFig2c | `h1_occ06_extfig2c_locale_ENVv2_cmm_pdb.tar.gz` (beads coloured by their v2 locale) | `code/chimerax/extfig2c_all.py` → `extfig2c_locale_scene.py` | V2 FINAL `ExtFig2c_renders` |

The cmm/pdb files were written on the cluster from the HSS files by `<Figure>/provenance/*.py` (not runnable here; needs the
HSS). Fig5b / ExtFig2c structures were regenerated with the ENV v2 locale table (`fig5_structures.py` /
`extfig2c_locale_structures.py` with `locale_seed20all_2p5r`). Fig5 tarballs hold SD + entropy only (no IQR).
`Fig4/run_chimerax.sh` was written for this package and has not been run (it opens a ChimeraX window); the Fig5 / ExtFig2
renderers are the package scripts (Fig5: one ChimeraX call per metric, SD / entropy only).

## Numbers that differ from the PDF (locale version of the PDF panel ≠ ENV v2)

| Panel | PDF | This package (ENV v2) | PDF was made with |
|---|---|---|---|
| Fig2d middle (ES chr2, entropy vs SD) | r = 0.47 | r = 0.429 (n 1,187) | ENV v1 (0.4718); old N565 = 0.393 |
| ExtFig2c bars, chr21:10.0–10.2 Mb box | entropy 0.95 | 0.90 (0.898); 0.69 / 2.16 unchanged | v1 (0.955) or old N565 (0.949) |
| ExtFig2d (ES, active / repressive) | r −0.14 (n 6,524) / +0.37 (n 6,791) | r **−0.1385 (n 6,528: spe 2,117, NL-act 4,411)** / **+0.3670 (n 6,787: env 4,323, NL-rep 2,464)** | ENV v1: −0.1362 / +0.3733 (reproduced exactly by `plot_extfig2d.py` with `EXTFIG2D_DATA=<dir with v1 tables>`) |
| ExtFig2e (ES chr11 / chr20) | r 0.035 / −0.593 | r −0.006 / −0.599 | ENV v1 (0.0346 / −0.5926); Fig5d of the PDF already shows the v2 chr11 value −0.006 |
| ExtFig6b (% spe-dominant in top 10 %) | 30.9 / 43.6 / 53.2 / 47.7 | 30.9 / 42.9 / 53.5 / 47.0 | old N565 locale (lamina ≤ 5 r) |

All other checked PDF numbers match ENV v2: Fig2c (−0.72 / 0.55 / −0.87), Fig2d box n (2,117 / 4,411 / 2,464 / 4,323),
Fig3e (−0.72 / −0.53 / −0.52 / −0.33), Fig5c, Fig5d, Fig5e (ES env 32.5 %), ExtFig2a (0.44 n 2,736 / 0.16 n 2,735),
ExtFig2b (−0.51; 0.37), ExtFig6a, ExtFig6c, ExtFig6d p values, ExtFig6e.

Layout notes: the PDF ExtFig6d x labels read "HCL" (should be HLC). Fig3e is drawn 1 x 4 (the PDF inserts an empty iHLC
column); Fig3c already contains the blank iHLC column.

**Text to fix in the manuscript:** the Fig. 2a legend says "200-kb resolution, 118 nm sphere radius". The ES occ0.6 model
uses an **ellipsoidal nucleus (semi-axes 4,368 x 3,369 x 2,889 nm)** and **bead radius 94.9 nm**, N = 29,838 beads
(2 x 14,385 autosomal + 781 chrX + 287 chrY).

## Panels not included

Experimental data or made by others / outside this package: Fig1 (all); Fig2a, 2f, 2g; Fig3a, 3b; the scHi-C rows of
Fig3c–e; Fig3f–h; Fig4c, 4d, 4g, 4h; Fig6 (all, incl. 6c / 6e, made outside this package); ExtFig1; ExtFig2f (contact maps,
made outside this package), ExtFig2g; ExtFig3; ExtFig4a–c, e–h and the plate / droplet rows of ExtFig4d; ExtFig5a;
ExtFig7 (all, incl. 7g, made outside this package).

## Verification (2026-11-12)

- `bash run_all.sh` from a clean output tree: all 8 figures exit 0 (≈45 s).
- Source data vs reference runs: 19 tables md5-identical; 11 tables = reference minus the IQR rows (every remaining line
  identical; ExtFig2a tables also have the droplet rows first); Fig5d / ExtFig2e / ExtFig6c Spearman tables = rows of `new_panels_spearman.tsv` (only the figure tag renamed);
  ExtFig6d box table = reference rows (order differs); ExtFig2b-right r identical; ExtFig2d reproduces the v1 table exactly
  when pointed at the v1 locale tables.
- PNGs: 31 figures pixel-identical to the reference PNGs. Changed layouts (Fig2b, Fig3c, ExtFig2a, ExtFig2b left, ExtFig2d,
  ExtFig6d) were checked by eye against the PDF.
- No absolute `/Users/yewang/plots` or `Figures_reproducible` path in any `code/`; no symlinks; `data/checksums.md5` per figure.
- Locale switch: after moving the tables to `data/locale/v2_ab4r/`, the default `run_all.sh` reproduced all 36 source-data
  tables byte-for-byte (md5). In a scratch copy, `set_locale_version.sh testv <v2 tables>` + `LOCALE_VERSION=testv
  bash run_all.sh` exited 0, wrote only `output_testv/` in the six locale folders, and gave the same 34 tables as the
  default; an unknown version stops with "no locale tables in data/locale/<version>".

- DISTCT_NORM switch (2026-11-12): with `SOURCE_DATE_EPOCH` fixed, the default `run_all.sh` after the change gave all 181
  output files md5-identical to the run before the change (PDFs differ otherwise only by their creation date);
  `DISTCT_NORM=r` and `=Req` runs exit 0 and write only `output_distCT-<norm>/`. Code backup before the change:
  `/Users/yewang/plots/distCT_normalisation/backup/Figures_final_code_before_DISTCT_NORM_2026-11-12.tar.gz`.

## Where the code came from

| Figure | Copied from | Edits |
|---|---|---|
| Fig2 | `~/Figures_reproducible/Fig2/code` (`plot_fig2d.py` from `~/plots/Figures_noHB_ENVv2/Fig2/code`, box order) | SD only; output names; IQR track removed from Fig2b |
| Fig3 | `~/Figures_reproducible/Fig3/code` (`plot_fig3d.py` from `Figures_noHB_ENVv2`, thick GMM lines) | SD only; IQR row / IQR & Mean histograms removed; names |
| Fig4 | `~/Figures_reproducible/Fig4/code` | SD only; Fig4k (ablation) not in the PDF, not copied; new `chimerax/fig4_render_sessions.py` |
| Fig5 | `~/Figures_reproducible/Fig5/code` (common, colour bars, ChimeraX); `~/plots/Figures_updated_ENV_locale/new_panels/code` (`plot_fig5ce_per_locale.py`; `plot_new_panels.py` split into `plot_fig5d.py`, `ExtFig2/code/plot_extfig2e.py`, `ExtFig6/code/plot_extfig6c.py`) | import path, names, IQR colour bar removed, ChimeraX runner loops per metric |
| ExtFig2 | `~/Figures_reproducible/ExtFig2/code`; `~/plots/sdrad_vs_sddistct_4stage/code`; new `plot_extfig2d.py` | SD only, 1 x 2 layouts, old 2e → 2h |
| ExtFig4 / 5 | `~/Figures_reproducible/ExtFig{4,5}/code` | SD only; names |
| ExtFig6 | `~/Figures_reproducible/ExtFig6/code` (common, 6ab, spe_topvar); `Figures_noHB_ENVv2/ExtFig6/code` (`plot_extfig6e_solid.py`; ExtFig6d drawing code) | SD only; names; new combined `plot_extfig6d.py` |

ES / DE were modelled from fmax16-normalised Hi-C, HB / HLC from fmax24; stage differences partly reflect modelling
parameters (see the earlier package READMEs).
