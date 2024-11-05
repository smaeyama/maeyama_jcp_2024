#!/bin/sh

### gamme_e scan for rotating flux-tube model (Lx=100)
cat ../../../f0.62_gamma0_salpha_nl/hst/gkvp.qes.0*     > qes_gamma0.00.dat
cat ../../../f0.62_gamma0.05_salpha_nl/hst/gkvp.qes.0.* > qes_gamma0.05.dat
cat ../../../f0.62_gamma0.1_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.10.dat
cat ../../../f0.62_gamma0.2_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.20.dat
cat ../../../f0.62_gamma0.3_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.30.dat
cat ../../../f0.62_gamma0.35_salpha_nl/hst/gkvp.qes.0.* > qes_gamma0.35.dat
cat ../../../f0.62_gamma0.4_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.40.dat
cat ../../../f0.62_gamma0.5_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.50.dat
cat ../../../f0.62_gamma0.6_salpha_nl/hst/gkvp.qes.0.*  > qes_gamma0.60.dat

### Lx scan (50<Lx<800) for no-shear-flow case (gamma_e=0)
cat ../../../f0.61_gamma0_salpha_nl_nx028mj2/hst/gkvp.qes.0.*  > qes_remap_gamma0.00_nx028mj02.dat
cat ../../../f0.61_gamma0_salpha_nl_nx042mj3/hst/gkvp.qes.0.*  > qes_remap_gamma0.00_nx042mj03.dat
cat ../../../f0.61_gamma0_salpha_nl/hst/gkvp.qes.0.*           > qes_remap_gamma0.00_nx055mj04.dat
cat ../../../f0.61_gamma0_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.00_nx110mj08.dat
cat ../../../f0.61_gamma0_salpha_nl_nx220mj16/hst/gkvp.qes.0.* > qes_remap_gamma0.00_nx220mj16.dat
cat ../../../f0.61_gamma0_salpha_nl_nx440mj32/hst/gkvp.qes.0.* > qes_remap_gamma0.00_nx440mj32.dat

### Lx scan (50<Lx<800) for rotating flux-tube model (gamma_e=0.05)
cat ../../../f0.62_gamma0.05_salpha_nl_nx028mj2/hst/gkvp.qes.0.*  > qes_gamma0.05_nx028mj02.dat
cat ../../../f0.62_gamma0.05_salpha_nl_nx042mj3/hst/gkvp.qes.0.*  > qes_gamma0.05_nx042mj03.dat
cat ../../../f0.62_gamma0.05_salpha_nl/hst/gkvp.qes.0.*           > qes_gamma0.05_nx055mj04.dat
cat ../../../f0.62_gamma0.05_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_gamma0.05_nx110mj08.dat
cat ../../../f0.62_gamma0.05_salpha_nl_nx220mj16/hst/gkvp.qes.0.* > qes_gamma0.05_nx220mj16.dat
cat ../../../f0.62_gamma0.05_salpha_nl_nx440mj32/hst/gkvp.qes.0.* > qes_gamma0.05_nx440mj32.dat

### Lx scan (50<Lx<800) for time-discontinuous wave-vector remap method (gamma_e=0.05)
cat ../../../f0.61_gamma0.05_salpha_nl_nx028mj2/hst/gkvp.qes.0.*  > qes_remap_gamma0.05_nx028mj02.dat
cat ../../../f0.61_gamma0.05_salpha_nl_nx042mj3/hst/gkvp.qes.0.*  > qes_remap_gamma0.05_nx042mj03.dat
cat ../../../f0.61_gamma0.05_salpha_nl/hst/gkvp.qes.0.*           > qes_remap_gamma0.05_nx055mj04.dat
cat ../../../f0.61_gamma0.05_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.05_nx110mj08.dat
cat ../../../f0.61_gamma0.05_salpha_nl_nx220mj16/hst/gkvp.qes.0.* > qes_remap_gamma0.05_nx220mj16.dat
cat ../../../f0.61_gamma0.05_salpha_nl_nx440mj32/hst/gkvp.qes.0.* > qes_remap_gamma0.05_nx440mj32.dat

### gamme_e scan for time-discontinuous wave-vector remap method (Lx=100)
cat ../../../f0.61_gamma0.1_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.10_nx055mj04.dat
cat ../../../f0.61_gamma0.2_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.20_nx055mj04.dat
cat ../../../f0.61_gamma0.3_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.30_nx055mj04.dat
cat ../../../f0.61_gamma0.35_salpha_nl/hst/gkvp.qes.0.* > qes_remap_gamma0.35_nx055mj04.dat
cat ../../../f0.61_gamma0.4_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.40_nx055mj04.dat
cat ../../../f0.61_gamma0.5_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.50_nx055mj04.dat
cat ../../../f0.61_gamma0.6_salpha_nl/hst/gkvp.qes.0.*  > qes_remap_gamma0.60_nx055mj04.dat

### gamme_e scan for time-discontinuous wave-vector remap method (Lx=200)
cat ../../../f0.61_gamma0.1_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.10_nx110mj08.dat
cat ../../../f0.61_gamma0.2_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.20_nx110mj08.dat
cat ../../../f0.61_gamma0.3_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.30_nx110mj08.dat
cat ../../../f0.61_gamma0.4_salpha_nl_nx110mj8/hst/gkvp.qes.0.*  > qes_remap_gamma0.40_nx110mj08.dat

### gamme_e scan for time-discontinuous wave-vector remap method (Lx=400)
cat ../../../f0.61_gamma0.1_salpha_nl_nx220mj16/hst/gkvp.qes.0.*  > qes_remap_gamma0.10_nx220mj16.dat
cat ../../../f0.61_gamma0.2_salpha_nl_nx220mj16/hst/gkvp.qes.0.*  > qes_remap_gamma0.20_nx220mj16.dat
cat ../../../f0.61_gamma0.3_salpha_nl_nx220mj16/hst/gkvp.qes.0.*  > qes_remap_gamma0.30_nx220mj16.dat
cat ../../../f0.61_gamma0.4_salpha_nl_nx220mj16/hst/gkvp.qes.0.*  > qes_remap_gamma0.40_nx220mj16.dat

### gamme_e scan for time-discontinuous wave-vector remap method (Lx=800)
cat ../../../f0.61_gamma0.1_salpha_nl_nx440mj32/hst/gkvp.qes.0.*  > qes_remap_gamma0.10_nx440mj32.dat
cat ../../../f0.61_gamma0.2_salpha_nl_nx440mj32/hst/gkvp.qes.0.*  > qes_remap_gamma0.20_nx440mj32.dat
cat ../../../f0.61_gamma0.3_salpha_nl_nx440mj32/hst/gkvp.qes.0.*  > qes_remap_gamma0.30_nx440mj32.dat
cat ../../../f0.61_gamma0.4_salpha_nl_nx440mj32/hst/gkvp.qes.0.*  > qes_remap_gamma0.40_nx440mj32.dat

