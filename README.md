# Model (IGM) panels of the single-cell genome-folding variability manuscript — figure code

> The plate- and droplet-based scHi-C cis/total annotation tables in `data/schic/` are the manuscript's Supplementary
> Tables 2 and 3; please cite the paper when using them.

This repository regenerates every **statistical plot of the IGM-model panels** of the main figures (Fig. 2–5) and Extended
Data figures (ExtFig. 2, 4, 5, 6), together with the **source-data table** behind each panel, from small pre-computed
per-bin tables. The 3D structure renders (ChimeraX) additionally need a Zenodo download (see below).

Models: IGM population models of four differentiation stages — ES, DE, HB, HLC (model ids `h1`, `DE`, `HB`, `mHEP`),
occupancy 0.6, 200 kb, 1000 diploid structures each. Defaults = what the manuscript uses: nuclear locale **ENV v2**
(`data/locale/v2_ab4r/`) and dist-to-CT normalisation **e2c** (signed distance / core centroid-to-surface distance,
`data/distCT/e2c/`). Definitions and per-panel notes: `docs/README_full.md`.

## Quick start (macOS or Linux, ~1 min after the environment is installed)

```bash
git clone https://github.com/alberlab/sc-genome-folding-variability-IGM-figures.git
cd sc-genome-folding-variability-IGM-figures
conda env create -f environment.yml          # env "sc-folding-figures" (conda-forge; python 3.13, pinned versions)
conda activate sc-folding-figures
bash run_all.sh                              # all 8 figure folders; or e.g. bash Fig2/run.sh for one figure
```

- Without conda: `python3.13 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt && bash run_all.sh`.
- Any interpreter with the packages of `environment.yml` can be used: `PYTHON=/path/to/python bash run_all.sh`.
- `run_all.sh` prints the Python and package versions it uses. Runtime: ~45 s on an Apple M-series Mac (conda env
  creation ~1–3 min).
- Exact environments the authors used (macOS arm64 `conda list --export` per figure): `envs/authors_osx-arm64/`.
- Input integrity: `md5sum -c data/checksums.md5` from the repo root (macOS without coreutils:
  `while read -r m f; do [ "$(md5 -q "$f")" = "$m" ] || echo "FAIL $f"; done < data/checksums.md5`).

## Outputs

`<Figure>/output/figures/<panel>_<description>.{png,pdf}` and `<Figure>/output/source_data/*.tsv` (not tracked by git;
created by the run). A fresh run reproduces the authors' 36 source-data tables and 45 PNGs byte-for-byte (PDFs differ only
in their embedded creation date).

## Layout

```
run_all.sh                 runs every <Figure>/run.sh
environment.yml            conda env (requirements.txt = pip equivalent)
data/locale/v2_ab4r/       nuclear-locale tables, 4 stages (shared)
data/distCT/e2c/           dist-to-CT per-bin stats (4 stages) + per-locus samples (.npz) + SD(RAD) vs SD(dist to CT) table
data/radial/               per-bin radial position, 4 stages
data/schic/                plate / droplet scHi-C cis/total annotation tables (unpublished; Supplementary Tables 2/3)
data/checksums.md5         md5 of every input table;  data/zenodo_manifest.tsv = files of the Zenodo download
<Figure>/code/             plot_<panel>.py per panel + <figure>_common.py; ChimeraX scripts in code/chimerax/
<Figure>/data/             input tables used only by this figure (Fig2, Fig4, ExtFig5)
<Figure>/run.sh            runs the panel scripts of this figure
<Figure>/provenance/       cluster scripts that computed the tables / structures from the IGM HSS files (not runnable here)
docs/README_full.md        detailed notes: definitions, reference checks, numbers that differ from the PDF, code origin
```

## Panel → script → output

Scripts are in `<Figure>/code/`, outputs in `<Figure>/output/figures/`.

| Panel | Content | Script | Output |
|---|---|---|---|
| Fig2b | chr11 1 Mb tracks: bulk Eig1, cis/total; plate / droplet Mean & SD(cisratio); model Mean & SD(dist to CT) | `plot_fig2b.py` | `Fig2b_chr11_tracks` |
| Fig2c | Mean vs Eig1 chr2 (colour SD) / SD vs Eig1 / SD vs Mean (colour Eig1), 1 Mb | `plot_fig2c.py` | `Fig2c_{left,middle,right}_*` |
| Fig2d | entropy by dominant locale (box order spe, NL-act, NL-rep, env) / entropy vs SD chr2 / entropy vs Eig1 200 kb | `plot_fig2d.py` | `Fig2d_{left,middle,right}_*` (+ `_alt_`: genome-wide SD, 1 Mb Eig1) |
| Fig2e | chr11 1 Mb tracks: locale entropy, RAD, SD(RAD) | `plot_fig2e.py` | `Fig2e_chr11_tracks` |
| Fig2h | model chr11 5 loci inside/outside bars + histogram | `plot_fig2h.py` | `Fig2h_chr11_5loci_bar_hist` |
| Fig3c | IGM rows: Mean & SD(dist to CT), chr2 + chr11, ES DE HB (iHLC blank) HLC | `plot_fig3c.py` | `Fig3c_IGM_all_stages_*` (+ one file per stage) |
| Fig3d | IGM columns: SD(dist to CT) histogram; entropy histogram + GMM (thick lines lw 8 / 6.5) | `plot_fig3d.py` | `Fig3d_IGM_SDdistCT_hist`, `Fig3d_IGM_entropy_hist_GMM` |
| Fig3e | IGM row: Mean(dist to CT) vs own-stage Eig1, chr2, colour SD (1 x 4; the PDF adds the blank iHLC column) | `plot_fig3e.py` | `Fig3e_IGM_Mean_vs_Eig1_chr2_colorbySD` |
| Fig4a / 4e / 4i | ChimeraX: nuclear shape half view / chr11 CT / Shannon index | `chimerax/fig4_render_sessions.py` (`run_chimerax.sh`) | `Fig4a/`, `Fig4e/`, `Fig4i/` (ChimeraX; Zenodo) |
| Fig4b | mean SD(dist to CT) per radial decile | `plot_fig4b.py` | `Fig4b_SDdistCT_vs_radial_decile` |
| Fig4f | chr11 CT volume ratio / compactness (packing density) violins | `plot_fig4f.py` | `Fig4f_chr11_{CT_volume_ratio,compactness}` |
| Fig4j | mean Shannon index per radial decile | `plot_fig4j.py` | `Fig4j_Shannon_vs_radial_decile` |
| Fig5a | ChimeraX half-sphere, SD(dist to CT) + speckle centroids | `chimerax/fig5ab_*.py` (`run_chimerax.sh SD`); colour bar `plot_fig5ab_colorbars.py` | `Fig5a/` (ChimeraX; Zenodo) |
| Fig5b | ChimeraX half-sphere, locale entropy (v2) | `chimerax/fig5ab_*.py` (`run_chimerax.sh entropy`) | `Fig5b/` (ChimeraX; Zenodo) |
| Fig5c | composition of all bead instances: spe, NL-act, nuc, NL-rep, env, other | `plot_fig5ce_per_locale.py` | `Fig5c_locale_composition_stacked` |
| Fig5d | entropy vs SD, chr2 / chr11 x 4 stages | `plot_fig5d.py` | `Fig5d_entropy_vs_SDdistCT_chr2_chr11` |
| Fig5e | 5-class dominant-locale Sankey ES → DE → HB → HLC | `plot_fig5ce_per_locale.py` | `Fig5e_dominant_locale_sankey` |
| ExtFig2a | model SD(dist to CT) vs droplet / plate SD(cisratio), 1 Mb density | `plot_extfig2a.py` | `ExtFig2a_model_vs_schic_SDcisratio` |
| ExtFig2b left | SD(dist to CT) vs RAD, ES density | `plot_extfig2b.py` | `ExtFig2b_left_SDdistCT_vs_RAD` |
| ExtFig2b right | SD(RAD) vs SD(dist to CT), ES density | `plot_extfig2b_right.py` | `ExtFig2b_right_SDRAD_vs_SDdistCT_ES` |
| ExtFig2c | ES locale renders (5 classes) + 3 example loci locale-frequency bars | `plot_extfig2c.py`; `chimerax/extfig2c_*.py` (`run_chimerax.sh`) | `ExtFig2c_example_loci_locale_freq`, `ExtFig2c/` (ChimeraX; Zenodo) |
| ExtFig2d | ES all autosomes, entropy vs SD: active (spe + NL-act dominant) / repressive (env + NL-rep + nuc), opacity = strength | `plot_extfig2d.py` (shared `data/`) | `ExtFig2d_ES_entropy_vs_SDdistCT_active_repressive` |
| ExtFig2e | ES chr11 / chr20 entropy vs SD | `plot_extfig2e.py` (shared `data/`) | `ExtFig2e_ES_entropy_vs_SDdistCT_chr11_chr20` |
| ExtFig2h | model chr2 3 loci inside/outside bars + histogram | `plot_extfig2h.py` | `ExtFig2h_chr2_3loci_bar_hist` |
| ExtFig4d | IGM row: chr11 5 loci, Mean & SD(dist to CT) across stages | `plot_extfig4.py` | `ExtFig4d_IGM_chr11_5loci_Mean_SDdistCT` |
| ExtFig5b | model chr2 CT volume ratio / compactness (packing density) violins | `plot_extfig5b.py` | `ExtFig5b_chr2_{CT_volume_ratio,compactness}` |
| ExtFig6a | top: SD(dist to CT) vs RAD; bottom: entropy vs RAD; 4 stages, colour = dominant locale | `plot_extfig6ab.py` | `ExtFig6a_{top_SDdistCT,bottom_entropy}_vs_RAD` |
| ExtFig6b | % of top-10 % SD(dist to CT) bins that are spe-dominant | `plot_extfig6_spe_topvar.py` | `ExtFig6b_top10_SDdistCT_fraction_spe_dominant` (+ `_alt_`: % of spe bins in top 10 %) |
| ExtFig6c | entropy vs SD, all autosomes / chr20, 4 stages | `plot_extfig6c.py` (shared `data/`) | `ExtFig6c_entropy_vs_SDdistCT_merged_chr20` |
| ExtFig6d | top: SD(dist to CT) boxes; bottom: entropy boxes + two-sided Mann–Whitney p (each stage vs ES); spe, NL-act, env, NL-rep | `plot_extfig6d.py` | `ExtFig6d_box_SD_entropy_by_dominant_locale_MWU` |
| ExtFig6e | entropy vs SD 4 x 4 (rows spe, env, NL-act, NL-rep), solid colours + legend | `plot_extfig6e_solid.py` | `ExtFig6e_entropy_vs_SD_by_dominant_locale_solid` |

Not included (experimental data or made outside this package): Fig1; Fig2a, 2f, 2g; Fig3a, 3b, the scHi-C rows of
Fig3c–e, Fig3f–h; Fig4c, 4d, 4g, 4h; Fig6; ExtFig1; ExtFig2f, 2g; ExtFig3; ExtFig4a–c, e–h and the scHi-C rows of 4d;
ExtFig5a; ExtFig7.

## ChimeraX renders (Fig4a/e/i, Fig5a/b, ExtFig2c)

Not produced by `run_all.sh`. They need (1) the Zenodo download — cmm/pdb structure tarballs, nuclear-envelope maps
(`*_nucleus_shell_pores.mrc`, up to 55 MB each) and the Fig4 ChimeraX sessions (`.cxs`) — fetched by `bash fetch_data.sh`
into the paths listed in `data/zenodo_manifest.tsv` (68 files, ~350 MB; **DOI to be added at publication**, the script is a
placeholder until then), and (2) **UCSF ChimeraX 1.11.1** (macOS arm64; the sessions were written by ChimeraX 1.11,
commit 60235d30). Then: `bash Fig4/run_chimerax.sh`, `bash Fig5/run_chimerax.sh [SD] [entropy]`,
`bash ExtFig2/run_chimerax.sh` (a ChimeraX window opens and closes by itself; the repo path must not contain spaces).

## Switches (optional, non-default inputs not shipped)

`LOCALE_VERSION=<v>` reads `data/locale/<v>/` and writes `<Figure>/output_<v>/`; `DISTCT_NORM=r|Req` reads
`data/distCT/<norm>/` and writes `<Figure>/output_distCT-<norm>/`. Install such tables with
`bash set_locale_version.sh <v> <src_dir>` / `bash set_distct_norm.sh <r|Req> <src_dir>` (formats in the script headers).
Only the manuscript defaults are in this repository; alternative tables are available from the authors on request.

## Data provenance

All model tables were computed on the UCLA Hoffman2 cluster from the IGM structure populations (`.hss`, not included;
scripts in `<Figure>/provenance/` and documented in `docs/README_full.md`): dist-to-CT statistics (CT core = convex hull of
the top 30 % strict-cis-cisratio beads per chromosome copy, contacts within 2.5 bead radii; 2 copies x 1000 structures per
200 kb bin), radial positions, nuclear-locale frequencies / entropy (speckle 5 r, nucleolus 6 r, envelope = outermost 20 %
of beads + 2.5 r, NL-act / NL-rep by A/B ratio within 4 r), CT volume / compactness and Shannon index. `Fig2/data/h1_bulk_*`
are bulk H1 Hi-C compartment (Eig1) and cis/total tracks; `data/schic/*.xlsx` are the collaborators' scHi-C annotation
tables (column E1 = bulk Eig1 of each stage). These tables will be deposited with the paper.
