#!/usr/bin/env python3
"""
Write the summary log of one isolate-wf sample: what every step did and what it wrote.

Each step's files are staged in a directory of their own by the SAMPLE_SUMMARY process; a
step that was skipped or never reached for this sample has an empty or missing directory,
and the log says which of the two it was. Only the standard library is used, so the script
runs in any container that has Python 3.
"""

import argparse
import csv
import re
from pathlib import Path

RULE = "-" * 80


# --------------------------------------------------------------------------- helpers


def files_in(directory, pattern="*"):
    """Files matching a pattern, searched recursively, in a stable order."""
    if directory is None or not Path(directory).is_dir():
        return []
    return sorted(p for p in Path(directory).rglob(pattern) if p.is_file())


def first_file(directory, pattern="*"):
    found = files_in(directory, pattern)
    return found[0] if found else None


def num(value, decimals=0):
    """A number with thousands separators; strings that are not numbers pass through."""
    try:
        value = float(str(value).replace(",", ""))
    except ValueError:
        return str(value)
    return f"{value:,.{decimals}f}"


def plural(count, word):
    return f"{count} {word}" + ("" if count == 1 else "s")


def pct(part, whole):
    return f"{100 * part / whole:.1f}%" if whole else "n/a"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_flat_yaml(path):
    """The top-level `key: value` pairs of a YAML file; nested blocks are ignored."""
    values = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^([A-Za-z_]\w*):\s*(\S.*)?$", line)
        if match and match.group(2) is not None:
            values[match.group(1)] = match.group(2).strip().strip('"')
    return values


def read_subsample_yaml(path):
    """subsample.yaml: the input read metrics, and a count/bases/n50 entry per subset."""
    values = read_flat_yaml(path)
    subsets = []
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^\s*(-\s+)?(count|bases|n50):\s*(\d+)", line)
        if not match:
            continue
        if match.group(1):
            subsets.append({})
        if subsets:
            subsets[-1][match.group(2)] = int(match.group(3))
    return values, subsets


def read_versions(path):
    """Tool versions from the collated versions file, as {tool: version}."""
    versions = {}
    if path is None or not Path(path).is_file():
        return versions
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^\s+([^:\s]+):\s*(.+)$", line)
        if match:
            versions[match.group(1).lower()] = match.group(2).strip().strip('"')
    return versions


def read_fasta_headers(path):
    """(name, length, header) for every record of a FASTA file."""
    records = []
    name, header, length = None, "", 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, length, header))
                header = line[1:]
                name = header.split()[0] if header.split() else ""
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        records.append((name, length, header))
    return records


def read_info(path):
    info = {}
    for line in Path(path).read_text().splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            info[key] = value
    return info


def as_list(value):
    return [item for item in (value or "").split(",") if item]


# --------------------------------------------------------------------------- the log


class Log:
    def __init__(self, versions):
        self.lines = []
        self.versions = versions
        self.step = 0

    def add(self, text="", indent=0):
        self.lines.append(" " * indent + text if text else "")

    def section(self, title, tool=None, summary=None):
        self.step += 1
        version = self.versions.get(tool.lower()) if tool else None
        label = f"{tool} {version}" if version else tool
        heading = f"[{self.step}] {title}" + (f" ({label})" if label else "")
        self.add()
        self.add(RULE)
        self.add(heading)
        self.add(RULE)
        if summary:
            self.add(summary, 4)

    def outputs(self, paths):
        paths = [p for p in paths if p]
        if not paths:
            return
        self.add("Outputs:", 4)
        for path in paths:
            self.add(path, 6)

    def text(self):
        return "\n".join(self.lines).rstrip() + "\n"


def not_run(log, reason):
    log.add(f"Not run: {reason}.", 4)


# --------------------------------------------------------------------------- sections


def section_input(log, info):
    log.section("Input reads")
    log.add(info.get("reads", "unknown"), 4)
    log.add(f"Read type: {info.get('read_type', 'unknown')}", 4)


def section_nanoplot(log, info, directory, stage, sample, prefix_suffix, subdir):
    log.section(f"Read quality {stage}", "NanoPlot")
    if info.get("skip_nanoplot"):
        not_run(log, f"skipped ({info['skip_nanoplot']})")
        return
    stats_file = first_file(directory, "*NanoStats.txt")
    if stats_file is None:
        not_run(log, "no NanoPlot statistics for this sample")
        return
    stats = {}
    for line in stats_file.read_text().splitlines():
        if ":" in line and not line.startswith(">"):
            key, value = line.split(":", 1)
            stats[key.strip()] = value.strip()
    log.add("Summarised the length and quality of the reads.", 4)
    log.add(
        f"Reads: {num(stats.get('Number of reads', '?'))}"
        f" | Bases: {num(stats.get('Total bases', '?'))}"
        f" | N50: {num(stats.get('Read length N50', '?'))} bp"
        f" | Mean length: {num(stats.get('Mean read length', '?'), 1)} bp"
        f" | Mean quality: Q{stats.get('Mean read quality', '?')}",
        4,
    )
    log.outputs([f"quality_assessment/nanoplot/{subdir}/{sample}{prefix_suffix}_NanoPlot-report.html"])


def section_porechop(log, info, directory, sample):
    log.section("Adapter trimming", "Porechop_ABI")
    if info.get("skip_porechop"):
        not_run(log, f"skipped ({info['skip_porechop']})")
        return
    log_file = first_file(directory, "*.log")
    if log_file is None:
        not_run(log, "no Porechop_ABI log for this sample")
        return
    text = re.sub(r"\x1b\[[0-9;]*m", "", log_file.read_text())

    log.add("Searched the reads for adapters, trimmed them from the read ends and", 4)
    log.add("discarded the reads with an adapter in their middle.", 4)

    # The adapters used follow this heading, one `name: sequence` per line
    adapters = []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() == "Trimming adapters from read ends":
            for adapter_line in lines[index + 1:]:
                match = re.match(r"^\s+(\S.*?):\s*([ACGTN]+)\s*$", adapter_line)
                if not match:
                    break
                adapters.append(match.group(1))
            break
    if adapters:
        inferred = [a for a in adapters if a.startswith("Adapter_")]
        known = [a for a in adapters if not a.startswith("Adapter_")]
        if inferred:
            log.add(f"Adapters inferred from the reads: {', '.join(inferred)}", 4)
        if known:
            log.add(f"Known adapters matched: {', '.join(known)}", 4)

    for end in ("start", "end"):
        match = re.search(rf"([\d,]+) / ([\d,]+) reads had adapters trimmed from their {end} \(([\d,]+) bp removed\)", text)
        if match:
            trimmed, total = int(match.group(1).replace(",", "")), int(match.group(2).replace(",", ""))
            log.add(f"Trimmed from the {end}: {num(trimmed)} / {num(total)} reads ({pct(trimmed, total)}), {match.group(3)} bp removed", 4)
    match = re.search(r"([\d,]+) / ([\d,]+) reads were discarded based on middle adapters", text)
    if match:
        log.add(f"Discarded for a middle adapter: {match.group(1)} / {match.group(2)} reads", 4)
    match = re.search(r"([\d,]+) / ([\d,]+) reads were split based on middle adapters", text)
    if match:
        log.add(f"Split at a middle adapter: {match.group(1)} / {match.group(2)} reads", 4)
    log.outputs([f"quality_control/porechop/{log_file.name}"])


def section_chopper(log, info, directory, sample):
    log.section("Length and quality filtering", "chopper")
    if info.get("skip_chopper"):
        not_run(log, f"skipped ({info['skip_chopper']}); the reads were assembled as given")
        return
    log_file = first_file(directory, "*.log")
    log.add(
        f"Kept the reads of at least {num(info.get('chopper_minlength', '?'))} bp"
        f" and mean quality Q{info.get('chopper_minq', '?')} or above.",
        4,
    )
    match = re.search(r"Kept (\d+) reads out of (\d+) reads", log_file.read_text()) if log_file else None
    if match:
        kept, total = int(match.group(1)), int(match.group(2))
        log.add(f"Kept {num(kept)} / {num(total)} reads ({pct(kept, total)}); {num(total - kept)} removed", 4)
    log.outputs([f"reads/post_qc/{sample}.choppered.fastq.gz"])


def section_genome_size(log, info):
    log.section("Genome size estimate", "Autocycler")
    size = info.get("genome_size")
    if not size:
        not_run(log, "no genome size estimate for this sample")
        return
    log.add(f"Estimated the genome size from the reads: {num(size)} bp.", 4)
    log.add("It sets the depth of the read subsets and guides the assemblers.", 4)


def section_subsample(log, directory):
    log.section("Read subsampling", "Autocycler")
    yaml_file = first_file(directory, "subsample.yaml")
    if yaml_file is None:
        not_run(log, "no subsampling metrics for this sample")
        return
    values, subsets = read_subsample_yaml(yaml_file)
    log.add(
        f"Split {num(values.get('input_read_count', '?'))} reads"
        f" ({num(values.get('input_read_bases', '?'))} bp, N50 {num(values.get('input_read_n50', '?'))} bp)"
        f" into {len(subsets)} subsets, one per assembly round:",
        4,
    )
    for index, subset in enumerate(subsets, 1):
        log.add(
            f"Subset {index:02d}: {num(subset.get('count', '?'))} reads,"
            f" {num(subset.get('bases', '?'))} bp, N50 {num(subset.get('n50', '?'))} bp",
            6,
        )


def section_assembly(log, info, directory):
    log.section("Assembly of every subset", "Autocycler")
    assemblers = as_list(info.get("assemblers"))
    subsets = as_list(info.get("subsamples"))
    skipped = as_list(info.get("skipped_assemblers"))
    log.add(f"Ran `autocycler helper` with each of: {', '.join(assemblers) or 'none'}.", 4)
    if skipped:
        log.add(f"Skipped with --skip_<assembler>: {', '.join(skipped)}.", 4)

    fastas = {path.stem: path for path in files_in(directory, "*.fasta")}
    succeeded = 0
    rows = []
    for assembler in assemblers:
        for subset in subsets:
            fasta = fastas.get(f"{assembler}_{subset}")
            records = read_fasta_headers(fasta) if fasta and fasta.stat().st_size else []
            if records:
                succeeded += 1
                circular = sum("circular=true" in header.lower() for _, _, header in records)
                status = f"{plural(len(records), 'contig')} ({circular} circular), {num(sum(l for _, l, _ in records))} bp"
            else:
                status = "FAILED, no contigs"
            rows.append(f"{assembler:<12} subset {subset:<4} {status}")
    total = len(assemblers) * len(subsets)
    log.add(f"{succeeded} / {total} assemblies produced contigs:", 4)
    for row in rows:
        log.add(row, 6)
    log.add("The assemblies are not published one by one: they are in the unitig graph", 4)
    log.add("of the archive (autocycler_out/input_assemblies.gfa).", 4)
    return succeeded


def section_clustering(log, info, directory, assemblies_ok):
    log.section("Compression and clustering", "Autocycler")
    yaml_file = first_file(directory, "clustering.yaml")
    if yaml_file is None:
        not_run(log, "no assembly produced contigs" if not assemblies_ok else "no clustering for this sample")
        return None
    values = read_flat_yaml(yaml_file)
    log.add("Compressed every assembly into one unitig graph, then grouped their contigs", 4)
    log.add("into putative replicons (clusters) and kept those that pass QC.", 4)
    passed, failed = int(values.get("pass_cluster_count", 0)), int(values.get("fail_cluster_count", 0))
    log.add(f"Clusters: {passed} passed QC, {failed} failed", 4)
    log.add(
        f"Contigs:  {values.get('pass_contig_count', '?')} in passing clusters,"
        f" {values.get('fail_contig_count', '?')} in failing ones",
        4,
    )
    score = values.get("overall_clustering_score")
    if score:
        log.add(f"Overall clustering score: {float(score):.3f} (1 is best)", 4)
    sample = info["sample"]
    log.outputs([
        "autocycler_clustering/clustering.tsv, clustering.newick, clustering.yaml",
        f"autocycler_clustering/{sample}_autocycler_out.tar.gz",
    ])
    return passed


def parse_list(value):
    """`[8, 3]` from `autocycler table` into ['8', '3']."""
    return [item.strip() for item in (value or "").strip("[]").split(",") if item.strip()]


def section_consensus(log, info, metrics_dir, consensus_dir, passed_clusters):
    log.section("Trimming, resolving and consensus", "Autocycler")
    sample = info["sample"]
    consensus = first_file(consensus_dir, "*.fasta")
    metrics_file = first_file(metrics_dir, "*_autocycler_metrics.tsv")
    if consensus is None and metrics_file is None:
        if passed_clusters is None:
            not_run(log, "no clustering for this sample")
        elif passed_clusters == 0:
            not_run(log, "no cluster passed QC; curate the clusters by hand from the archive")
        else:
            not_run(log, "no consensus for this sample")
        return []

    log.add("Trimmed the overlaps of each passing cluster, resolved its repeats, and", 4)
    log.add("combined the clusters into the consensus assembly.", 4)
    metrics = read_tsv(metrics_file)[0] if metrics_file else {}
    sizes = parse_list(metrics.get("trimmed_cluster_size"))
    medians = parse_list(metrics.get("trimmed_cluster_median"))
    for index, size in enumerate(sizes, 1):
        median = medians[index - 1] if index - 1 < len(medians) else "?"
        log.add(f"Cluster {index}: {size} sequences after trimming, median length {num(median)} bp", 6)
    if metrics:
        resolved = metrics.get("consensus_assembly_fully_resolved", "?").lower() == "true"
        log.add(
            f"Consensus: {num(metrics.get('consensus_assembly_bases', '?'))} bp in"
            f" {metrics.get('consensus_assembly_unitigs', '?')} unitigs,"
            f" {'fully resolved' if resolved else 'NOT fully resolved (inspect consensus_assembly.gfa)'}",
            4,
        )

    records = read_fasta_headers(consensus) if consensus else []
    if records:
        log.add(f"Contigs of consensus_assembly.fasta ({len(records)}):", 4)
        for name, length, header in records:
            topology = "circular" if "circular=true" in header.lower() else "linear"
            log.add(f"{name:<10} {num(length):>12} bp  {topology}", 6)
    log.outputs([
        "autocycler_final/consensus_assembly.fasta" if consensus else None,
        "autocycler_final/consensus_assembly.gfa" if consensus else None,
        f"autocycler_final/{metrics_file.name}" if metrics_file else None,
        "../autocycler_metrics.tsv (every sample)" if metrics_file else None,
    ])
    return records


def section_dnaapler(log, info, directory, has_consensus):
    log.section("Reorientation", "Dnaapler")
    if info.get("skip_dnaapler"):
        not_run(log, f"skipped ({info['skip_dnaapler']})")
        return
    if not has_consensus:
        not_run(log, "no consensus assembly")
        return
    summary = first_file(directory, "*_reorientation_summary.tsv")
    log.add("Rotated the circular contigs to start at dnaA (chromosome), repA (plasmid)", 4)
    log.add("or terL (phage); linear contigs were left as they are.", 4)
    if summary:
        for row in read_tsv(summary):
            name = (row.get("Contig") or "?").split()[0]
            gene = row.get("Gene_Reoriented", "")
            if gene == "Contig_ignored":
                outcome = "linear, left as is"
            elif gene and gene.lower() not in ("none", "no_hit", "nan"):
                identity = row.get("Identity_Percentage", "")
                outcome = f"now starts at {gene}" + (f" ({row.get('Top_Hit', '')}, {identity}% identity)" if identity else "")
            else:
                outcome = "no dnaA, repA or terL found, left as is"
            log.add(f"{name:<10} {outcome}", 6)
    log.outputs([f"dnaapler/dnaapler_out/{info['sample']}_reoriented.fasta"])


def section_checkm2(log, info, directory, has_consensus):
    log.section("Completeness and contamination", "CheckM2")
    if info.get("skip_checkm2"):
        not_run(log, f"skipped ({info['skip_checkm2']})")
        return
    if not has_consensus:
        not_run(log, "no consensus assembly")
        return
    report = first_file(directory, "quality_report.tsv")
    if report is None:
        not_run(log, "no CheckM2 report for this sample")
        return
    rows = read_tsv(report)
    log.add("Estimated the completeness and contamination of the consensus assembly", 4)
    log.add("from its marker genes.", 4)
    for row in rows:
        log.add(
            f"Completeness: {row.get('Completeness', '?')}% | Contamination: {row.get('Contamination', '?')}%"
            f" | Model: {row.get('Completeness_Model_Used', '?')}",
            4,
        )
    log.outputs(["checkm2/checkm2_out/quality_report.tsv"])


def section_genomad(log, info, directory, has_consensus):
    log.section("Plasmids, viruses and proviruses", "geNomad")
    if info.get("skip_genomad"):
        not_run(log, f"skipped ({info['skip_genomad']})")
        return
    if not has_consensus:
        not_run(log, "no consensus assembly")
        return
    plasmid_file = first_file(directory, "*_plasmid_summary.tsv")
    virus_file = first_file(directory, "*_virus_summary.tsv")
    if plasmid_file is None and virus_file is None:
        not_run(log, "no geNomad summary for this sample")
        return
    log.add("Classified the contigs, and regions of contigs, as plasmid or virus.", 4)
    plasmids = read_tsv(plasmid_file) if plasmid_file else []
    viruses = read_tsv(virus_file) if virus_file else []
    proviruses = [row for row in viruses if row.get("topology", "").lower() == "provirus"]
    viruses = [row for row in viruses if row not in proviruses]
    for label, rows in (("Plasmids", plasmids), ("Viruses", viruses), ("Proviruses", proviruses)):
        log.add(f"{label}: {len(rows)}", 4)
        for row in rows:
            where = f", {row['coordinates']}" if row.get("coordinates") and row["coordinates"] != "NA" else ""
            score = row.get("plasmid_score") or row.get("virus_score") or "?"
            log.add(
                f"{row.get('seq_name', '?')}: {num(row.get('length', '?'))} bp,"
                f" {row.get('topology', '?')}{where}, score {score}",
                6,
            )
    prefix = plasmid_file.parent.name if plasmid_file else virus_file.parent.name
    log.outputs([f"genomad/{prefix}/"])


# --------------------------------------------------------------------------- main


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--info", required=True, help="TSV of key/value pairs from the workflow")
    parser.add_argument("--versions", help="Collated software versions YAML")
    parser.add_argument("--inputs", required=True, help="Directory holding one subdirectory per step")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    info = read_info(args.info)
    sample = info["sample"]
    inputs = Path(args.inputs)
    log = Log(read_versions(args.versions))

    title = f"isolate-wf {info.get('pipeline_version', '')} - summary of sample {sample}".replace("  ", " ")
    log.add(title)
    log.add("=" * len(title))
    log.add(f"Paths are relative to sample_{sample}/, the folder of this log.")
    result_index = len(log.lines)

    section_input(log, info)
    section_nanoplot(log, info, inputs / "nanoplot_raw", "before QC", sample, "_raw", "raw")
    section_porechop(log, info, inputs / "porechop", sample)
    section_chopper(log, info, inputs / "chopper", sample)
    section_nanoplot(log, info, inputs / "nanoplot_qc", "after QC", sample, "_qc", "post_qc")
    section_genome_size(log, info)
    section_subsample(log, inputs / "subsample")
    assemblies_ok = section_assembly(log, info, inputs / "assemblies")
    passed = section_clustering(log, info, inputs / "clustering", assemblies_ok)
    records = section_consensus(log, info, inputs / "metrics", inputs / "consensus", passed)
    section_dnaapler(log, info, inputs / "dnaapler", bool(records))
    section_checkm2(log, info, inputs / "checkm2", bool(records))
    section_genomad(log, info, inputs / "genomad", bool(records))

    # The outcome goes at the top, where it is read first
    if records:
        circular = sum("circular=true" in header.lower() for _, _, header in records)
        outcome = (
            f"Consensus assembly: {plural(len(records), 'contig')} ({circular} circular),"
            f" {num(sum(length for _, length, _ in records))} bp"
        )
    elif not assemblies_ok:
        outcome = "No consensus assembly: every assembler failed."
    elif passed == 0:
        outcome = "No consensus assembly: no cluster passed Autocycler QC."
    else:
        outcome = "No consensus assembly."
    log.lines[result_index:result_index] = ["", f"Result: {outcome}"]

    if log.versions:
        log.add()
        log.add(RULE)
        log.add("Software versions of this run")
        log.add(RULE)
        for tool, version in sorted(log.versions.items()):
            log.add(f"{tool}: {version}", 4)

    Path(args.output).write_text(log.text())


if __name__ == "__main__":
    main()
