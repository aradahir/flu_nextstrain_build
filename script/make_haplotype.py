#!/usr/bin/env python3
import argparse
from collections import Counter
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--metadata", required=True)
p.add_argument("--output", required=True)
p.add_argument("--clade-column", default="clade")
p.add_argument("--subs-column", default="aaSubstitutions")
p.add_argument("--gene-prefix", default="HA1")
p.add_argument("--min-count", type=int, default=5)
p.add_argument("--founder-frac", type=float, default=0.9)
args = p.parse_args()

df = pd.read_csv(args.metadata, sep="\t", dtype=str)
clade = df[args.clade_column].fillna("unassigned")

prefix = args.gene_prefix + ":"
muts = df[args.subs_column].fillna("").apply(
    lambda s: {m.strip() for m in s.split(",") if m.strip().startswith(prefix)}
)

# clade-defining mutations = present in >= founder_frac of the clade's members
founder = {}
for c in clade.unique():
    sets = muts[clade == c]
    cnt = Counter(m for s in sets for m in s)
    founder[c] = {m for m, n in cnt.items() if n / len(sets) >= args.founder_frac}

def label(c, s):
    derived = sorted(s - founder[c], key=lambda m: int("".join(ch for ch in m.split(":")[1] if ch.isdigit()) or 0))
    return c if not derived else c + ":" + ",".join(m.split(":", 1)[1] for m in derived)

hap = pd.Series([label(c, s) for c, s in zip(clade, muts)], index=df.index)

# collapse rare haplotypes so the Auspice legend stays usable
counts = hap.value_counts()
rare = ~hap.isin(counts[counts >= args.min_count].index)
hap[rare] = clade[rare] + ":other"

df["haplotype"] = hap
df.to_csv(args.output, sep="\t", index=False)
print(f"{len(df)} rows, {df['haplotype'].nunique()} haplotypes")