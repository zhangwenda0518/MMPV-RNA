#!/usr/bin/env python3
"""Build plant virus full info with taxonomy + cluster + abundance from COBRA coverage."""
import argparse, csv, os
from collections import defaultdict
from pathlib import Path


def load_taxonomy(path):
    tax = {}
    if not Path(path).is_file(): return tax
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cid = row.get("contig_id", "").strip().strip('"')
            tax[cid] = {k: row.get(k, "") for k in
                ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species",
                 "confidence","primary_tool","Nucleic_acid"]}
    return tax


def load_clusters(path):
    """加载 vclust_clusters.tsv → {contig_id: cluster_rep_id} + 统计每簇大小"""
    cmap = {}
    cluster_sizes = defaultdict(int)
    if not Path(path).is_file(): return cmap, cluster_sizes
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                cid, rep = parts[0], parts[1]
                cmap[cid] = rep
                cluster_sizes[rep] += 1
    return cmap, cluster_sizes


def load_coverage(cobra_dir):
    """加载 COBRA coverage → {contig_id: {sample: coverage}}"""
    cov = defaultdict(dict)
    for cov_file in Path(cobra_dir).glob("*_clean/cobra_penguin_result/*.coverage.txt"):
        sample = cov_file.parts[-3].replace("_clean", "")
        for line in open(cov_file):
            parts = line.strip().split()
            if len(parts) >= 2:
                cov[parts[0]][sample] = float(parts[1])
    return cov


def load_ids(fasta_path):
    ids = set()
    if not Path(fasta_path).is_file(): return ids
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def coverage_stats(cs):
    """cs = {sample: coverage}. 返回 (n, max, min, mean)"""
    if not cs: return 0, 0, 0, 0
    vals = list(cs.values())
    n = len(vals)
    return n, max(vals), min(vals), sum(vals) / n


def is_reference_contig(cid):
    """检测是否为参考序列 (NCBI accession 格式, 如 OR489165.1, NC_001234.1)"""
    import re
    return bool(re.match(r'^[A-Z]{2,4}_\d+\.\d+$', cid)) or cid.startswith(('NC_', 'ref|'))


def main():
    p = argparse.ArgumentParser(description="Build plant virus full info")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--plant-fasta", required=True)
    p.add_argument("--cobra-dir", required=True)
    p.add_argument("--clusters", default=None, help="vclust_clusters.tsv")
    p.add_argument("--rescue-fa", default=None)
    args = p.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    tax = load_taxonomy(args.taxonomy)
    cov = load_coverage(args.cobra_dir)
    plant_ids = load_ids(args.plant_fasta)
    rescued_ids = load_ids(args.rescue_fa) if args.rescue_fa else set()
    cmap, cluster_sizes = load_clusters(args.clusters) if args.clusters else ({}, {})
    all_samples = sorted(set(s for cs in cov.values() for s in cs))

    n_plant = len(plant_ids)
    n_rescued = len(rescued_ids & plant_ids)

    # Per-contig TSV
    out_tsv = out / "All_plant.viruses_info.tsv"
    cols = ["contig_id", "cluster_id", "cluster_size", "source",
            "Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species",
            "confidence","primary_tool","rescued",
            "n_samples","max_cov","min_cov","mean_cov","ref_cov",
            "Nucleic_acid"]
    with open(out_tsv, "w") as tf:
        tf.write("\t".join(cols) + "\n")
        for cid in sorted(plant_ids):
            t = tax.get(cid, {})
            cs = cov.get(cid, {})
            rep = cmap.get(cid, cid)
            csize = cluster_sizes.get(rep, 1)
            n_samp, mx_cov, mn_cov, mn_cov_val = coverage_stats(cs)
            # ref coverage: coverage of cluster representative
            ref_cs = cov.get(rep, {})
            _, _, _, ref_cov = coverage_stats(ref_cs)
            is_rescued = "Y" if cid in rescued_ids else "N"
            source = "reference" if is_reference_contig(cid) else "contig"

            vals = [cid, rep, str(csize), source,
                    t.get("Realm",""), t.get("Kingdom",""), t.get("Phylum",""),
                    t.get("Class",""), t.get("Order",""), t.get("Family",""),
                    t.get("Genus",""), t.get("Species",""),
                    t.get("confidence",""), t.get("primary_tool",""), is_rescued,
                    str(n_samp), "{:.2f}".format(mx_cov), "{:.2f}".format(mn_cov),
                    "{:.2f}".format(mn_cov_val), "{:.2f}".format(ref_cov),
                    t.get("Nucleic_acid","")]
            tf.write("\t".join(vals) + "\n")

    # Genus summary
    genus_stats = defaultdict(lambda: {"n_contigs":0, "n_rescued":0, "sum_cov":0.0})
    for cid in plant_ids:
        g = tax.get(cid, {}).get("Genus", "Unclassified")
        genus_stats[g]["n_contigs"] += 1
        if cid in rescued_ids: genus_stats[g]["n_rescued"] += 1
        cs = cov.get(cid, {})
        if cs: genus_stats[g]["sum_cov"] += sum(cs.values()) / len(cs)
    with open(out / "all_plant_viruses_genus_summary.tsv", "w") as gf:
        gf.write("Genus\tn_contigs\tn_rescued\tavg_cov\n")
        for g, info in sorted(genus_stats.items(), key=lambda x: -x[1]["n_contigs"]):
            avg = info["sum_cov"] / info["n_contigs"] if info["n_contigs"] else 0
            gf.write("{}\t{}\t{}\t{:.2f}\n".format(g, info["n_contigs"], info["n_rescued"], avg))

    # Cluster summary
    cluster_out = out / "plant_virus_cluster_summary.tsv"
    with open(cluster_out, "w") as cf:
        cf.write("cluster_id\tn_contigs_plant\tn_members_total\tref_cov\tmax_cov_in_cluster\n")
        plant_clusters = defaultdict(set)
        for cid in plant_ids:
            rep = cmap.get(cid, cid)
            plant_clusters[rep].add(cid)
        for rep, members in sorted(plant_clusters.items(), key=lambda x: -len(x[1])):
            csize = cluster_sizes.get(rep, 1)
            ref_cs = cov.get(rep, {})
            _, _, _, ref_cv = coverage_stats(ref_cs)
            # max coverage of any plant member in this cluster
            max_in_cluster = max(
                (max(cov.get(m, {}).values()) if cov.get(m, {}) else 0)
                for m in members
            )
            cf.write("{}\t{}\t{}\t{:.2f}\t{:.2f}\n".format(
                rep, len(members), csize, ref_cv, max_in_cluster))

    print("Plant virus contigs: {}".format(n_plant))
    print("Rescued HQ: {}".format(n_rescued))
    print("Genera: {}".format(len(genus_stats)))
    print("Clusters: {}".format(len(plant_clusters)))
    print("Samples: {}".format(len(all_samples)))
    print("Output: {}".format(out_tsv))


if __name__ == "__main__":
    main()
