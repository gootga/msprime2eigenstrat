# build_topology_to_eigenstrat.py
import msprime
import numpy as np
import matplotlib.pyplot as plt
import tskit

# -------------------------
# Parameters
# -------------------------
n_pops = 15
inds_per_pop = 10        # diploid individuals per population
ploidy = 2
sequence_length = 10_000_000
recombination_rate = 1e-8
mutation_rate = 1.5e-8
seed_anc = 42
seed_mut = 43
out_prefix = "Demografia15pop"      # prefix for EIGENSTRAT output
# -------------------------

# -------------------------
# 1) Build custom Demography
# -------------------------
t_root = 2000    # ANC -> L + R
t_L_split = 1600 # L -> L1 + L2
t_L1 = 1400      # L1 -> POP1 + POP2
t_L2a = 1250     # POP3 vs rest
t_L2b = 1200     # POP4 vs rest
t_L2c = 1195     # POP5 vs rest
t_L2d = 1190     # POP6 vs POP7
t_R_split = 1500 # R -> R1 + R2
t_R1 = 1350      # R1 -> POP8 + POP9
t_R2a = 1100     # POP10 vs rest
t_R2b = 1080     # POP11 vs rest
t_R2c = 1060     # POP12 vs rest
t_R2d = 1040     # POP13 vs POP14

dem = msprime.Demography()
dem.add_population(name="ANC", initial_size=10000)

# Internal nodes
for internal in ["L", "R", "L1", "L2", "R1", "R2",
                 "L2a", "L2b", "L2c", "R2a", "R2b", "R2c"]:
    dem.add_population(name=internal, initial_size=10000)

# Terminal tips
for i in range(1, 16):
    dem.add_population(name=f"POP{i}", initial_size=10000)

# Splits
dem.add_population_split(time=t_root, derived=["L", "R"], ancestral="ANC")
dem.add_population_split(time=t_L_split, derived=["L1", "L2"], ancestral="L")
dem.add_population_split(time=t_L1, derived=["POP1", "POP2"], ancestral="L1")
dem.add_population_split(time=t_L2a, derived=["POP3", "L2a"], ancestral="L2")
dem.add_population_split(time=t_L2b, derived=["POP4", "L2b"], ancestral="L2a")
dem.add_population_split(time=t_L2c, derived=["POP5", "L2c"], ancestral="L2b")
dem.add_population_split(time=t_L2d, derived=["POP6", "POP7"], ancestral="L2c")
dem.add_population_split(time=t_R_split, derived=["R1", "R2"], ancestral="R")
dem.add_population_split(time=t_R1, derived=["POP8", "POP9"], ancestral="R1")
dem.add_population_split(time=t_R2a, derived=["POP10", "R2a"], ancestral="R2")
dem.add_population_split(time=t_R2b, derived=["POP11", "R2b"], ancestral="R2a")
dem.add_population_split(time=t_R2c, derived=["POP12", "R2c"], ancestral="R2b")
dem.add_population_split(time=t_R2d, derived=["POP13", "POP14"], ancestral="R2c")

# Admixture: POP15 from POP7 (20%) + POP8 (80%)
dem.add_admixture(time=50, derived="POP15", ancestral=["POP7", "POP8"], proportions=[0.2, 0.8])

dem.sort_events()

# -------------------------
# 2) Quick schematic tree (SVG + TXT)
# -------------------------
plot_samples = [msprime.SampleSet(1, population=f"POP{i}") for i in range(1, 16)]
ts_plot = msprime.sim_ancestry(
    samples=plot_samples,
    demography=dem,
    sequence_length=1,
    recombination_rate=0,
    random_seed=1
)

tree = ts_plot.first()
svg = tree.draw(format="svg")
with open("topology_schematic.svg", "w") as f:
    f.write(svg)

print("Saved schematic tree as topology_schematic.svg")

newick = tree.newick()
print("Tree (Newick):\n", newick)

with open("topology_schematic.txt", "w") as f:
    f.write(newick + "\n")

print("Saved Newick tree as topology_schematic.txt")

# -------------------------
# 3) Full ancestry + mutations
# -------------------------
samples = [
    msprime.SampleSet(num_samples=inds_per_pop, population=f"POP{i}", ploidy=ploidy)
    for i in range(1, n_pops + 1)
]

ts = msprime.sim_ancestry(
    samples=samples,
    demography=dem,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    random_seed=seed_anc,
)

mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed_mut, model="jc69")

# -------------------------
# 4) Build individuals & names
# -------------------------
hap_nodes = list(mts.samples())  # haploid samples in order
n_individuals = n_pops * inds_per_pop
assert len(hap_nodes) == n_individuals * ploidy

ind_node_groups = [hap_nodes[i:i+ploidy] for i in range(0, len(hap_nodes), ploidy)]
sample_names = []
for pop in range(1, n_pops + 1):
    for j in range(inds_per_pop):
        sample_names.append(f"POP{pop}_ind{j+1}")

# -------------------------
# 5) Extract variants & write EIGENSTRAT
# -------------------------
G_hap = mts.genotype_matrix()
variants = list(mts.variants())

geno_lines = []
snp_lines = []
for vi, var in enumerate(variants):
    hap_gt = G_hap[vi]
    derived_counts = np.array([hap_gt[i] + hap_gt[i+1] for i in range(0, len(hap_gt), 2)])
    geno_chars = []
    for dc in derived_counts:
        if int(dc) < 0:            # missing data
            geno_chars.append("9")
        elif dc in (0, 1, 2):      # valid diploid genotype
            ref = ploidy - int(dc) # convert derived count → reference count
            geno_chars.append(str(ref))
        else:
            # Multi-allelic or unexpected genotype → mark as missing
            geno_chars.append("9")
    geno_lines.append("".join(geno_chars))

    alleles = tuple(var.alleles)
    ref = alleles[0] if len(alleles) >= 1 else "N"
    alt = alleles[1] if len(alleles) >= 2 else "N"
    snp_id = f"rs{vi+1}"
    chrom = "1"
    gen_dist = "0"
    pos = int(var.site.position)
    snp_lines.append(f"{snp_id}\t{chrom}\t{gen_dist}\t{pos}\t{ref}\t{alt}")

with open(out_prefix + ".geno", "w") as f:
    f.write("\n".join(geno_lines) + "\n")

with open(out_prefix + ".snp", "w") as f:
    f.write("\n".join(snp_lines) + "\n")

with open(out_prefix + ".ind", "w") as f:
    for name in sample_names:
        pop_label = name.split("_")[0]
        f.write(f"{name}\tU\t{pop_label}\n")

print(f"Wrote {out_prefix}.geno  {out_prefix}.snp  {out_prefix}.ind  (SNPs={len(geno_lines)})")

