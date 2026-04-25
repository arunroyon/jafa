# Demo Commands

Create an 11.2 um PAH map:

```bash
jafa make-map \
  --cubes cube_ch1.fits cube_ch2.fits cube_ch3.fits \
  --feature PAH_11p2 \
  --outdir results/pah_11p2
```

Create a PAH ratio map:

```bash
jafa make-ratio \
  --cubes cube_ch1.fits cube_ch2.fits cube_ch3.fits \
  --feature1 PAH_11p2 \
  --feature2 PAH_6p2 \
  --snr-threshold 3 \
  --outdir results/ratio_11p2_6p2
```

Create a custom C60 18.9 um map:

```bash
jafa make-map \
  --cubes cube_ch3.fits \
  --wavelength 18.9 \
  --feature-window 18.82 19.02 \
  --anchors 18.55 18.70 19.15 19.30 \
  --outdir results/c60_18p9
```
