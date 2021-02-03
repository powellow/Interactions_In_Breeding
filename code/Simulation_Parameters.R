### Reference Population of Genotypes (RPG) & Genome ----

n_founders=1000 #Change this value to show effect of sampling on covariance/LD between QTL
n_qtl=1
n_chr=10
start_allele_freq = 0.8
n_traits = 3

### QTL Effects ----
add_qtl_effect_focal =
add_qtl_effect_background =
dom_qtl_effect_focal =
dom_qtl_effect_background =
epi_qtl_effect_focal =
epi_qtl_effect_background =

### Simulation Parameters ----
n_reps = 5
trait = 1
n_cycles = 20
nIndSel = 100
nCrosses = 100
nProgeny = 10
selection = "genotypic_values"

