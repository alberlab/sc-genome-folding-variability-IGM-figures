#!/usr/bin/env python3
"""Figure 2e: chr11 model tracks at 1 Mb — nuclear locale entropy, radial position (RAD) and SD(radial).
RAD = sqrt((x/a)^2 + (y/b)^2 + (z/c)^2) with the model's ellipsoidal nuclear semi-axes (a, b, c = 4368, 3369, 2889 nm),
averaged over the two copies and 1000 structures (0 = centre, 1 = envelope). SD(radial) = SD of RAD over the 2 copies x
1000 structures (2000 values per 200 kb bin, as for SD(dist to CT)). Same 1 Mb grid and mask as Fig. 2b.
Outputs: output/figures/Fig2e_chr11_tracks.*, output/source_data/Fig2e_chr11_1Mb.tsv
"""
from fig2_common import chr_tracks_1mb, plot_tracks, save_table

CHR = "chr11"
t, spans = chr_tracks_1mb(CHR)
save_table(t.reset_index().rename(columns={"index": "start"}).assign(chr=CHR)[
    ["chr", "start", "end", "masked", "entropy_bit", "norm_rad", "sd_rad"]], "Fig2e_chr11_1Mb.tsv")
TRACKS = [("entropy_bit", "Model nuclear locale entropy", "Entropy (bits)", "black", None),
          ("norm_rad", "Model radial position (RAD)", "RAD", "black", None),
          ("sd_rad", "Model SD(radial position)", "SD(RAD)", "red", None)]
plot_tracks(t, spans, TRACKS, CHR, f"H1 model tracks for {CHR} (1 Mb)", f"Fig2e_{CHR}_tracks")
