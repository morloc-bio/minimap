# Simulated genome and reads, and the minimap2 2.30 command line as the
# oracle for the minimap suite.
import os
import shutil
import random
import subprocess
import tempfile

_dir = None


def tmp_path(name):
    global _dir
    if _dir is None:
        _dir = tempfile.mkdtemp(prefix="minimap-test-")
    return os.path.join(_dir, name)


def _rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _mutate(rng, s, rate):
    out = []
    for c in s:
        u = rng.random()
        if u < rate:
            out.append(rng.choice([x for x in "ACGT" if x != c]))
        elif u < rate * 1.1:
            continue
        elif u < rate * 1.2:
            out.append(c + rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def write_inputs(seed):
    """A genome of three contigs (one holding a two-exon gene with a GT-AG
    intron), single reads from the genome (some crossing the intron as
    spliced reads) and read pairs. Returns the directory."""
    rng = random.Random(seed)
    d = tmp_path("in%d" % seed)
    os.makedirs(d, exist_ok=True)
    contigs = {"chrA": 60000, "chrB": 40000, "chrC": 20000}
    seqs = {n: "".join(rng.choice("ACGT") for _ in range(l)) for n, l in contigs.items()}
    a = seqs["chrA"]
    # exon 1 [20000, 20300), intron [20300, 22300), exon 2 [22300, 22600)
    a = a[:20300] + "GT" + a[20302:22298] + "AG" + a[22300:]
    seqs["chrA"] = a
    with open(os.path.join(d, "genome.fa"), "w") as f:
        for n, s in seqs.items():
            f.write(">%s\n" % n)
            for i in range(0, len(s), 60):
                f.write(s[i : i + 60] + "\n")
    names = list(seqs)
    reads = []
    for i in range(300):
        n = rng.choice(names)
        s = seqs[n]
        p = rng.randrange(0, len(s) - 150)
        r = _mutate(rng, s[p : p + 150], 0.01)
        if rng.random() < 0.5:
            r = _rc(r)
        reads.append(("u%d" % i, r))
    tx = a[20000:20300] + a[22300:22600]
    for i in range(100):
        p = rng.randrange(150, 300 - 100)
        r = tx[p : p + 150] if p + 150 <= len(tx) else tx[p:]
        r = _mutate(rng, r, 0.005)
        if rng.random() < 0.5:
            r = _rc(r)
        reads.append(("s%d" % i, r))
    with open(os.path.join(d, "reads.fq"), "w") as f:
        for h, r in reads:
            f.write("@%s\n%s\n+\n%s\n" % (h, r, "I" * len(r)))
    with open(os.path.join(d, "r_1.fq"), "w") as f1, open(os.path.join(d, "r_2.fq"), "w") as f2:
        for i in range(200):
            n = rng.choice(names)
            s = seqs[n]
            ins = rng.randint(250, 450)
            p = rng.randrange(0, len(s) - ins)
            frag = s[p : p + ins]
            if rng.random() < 0.5:
                frag = _rc(frag)
            m1, m2 = _mutate(rng, frag[:100], 0.01), _mutate(rng, _rc(frag)[:100], 0.01)
            f1.write("@p%d\n%s\n+\n%s\n" % (i, m1, "I" * len(m1)))
            f2.write("@p%d\n%s\n+\n%s\n" % (i, m2, "I" * len(m2)))
    return d


def _minimap2():
    p = shutil.which("minimap2")
    if p is None:
        raise RuntimeError("no minimap2 on the PATH (minimap/test/package.yaml tool-deps)")
    return p


def _paf(text):
    """PAF records as comparable tuples: the 12 columns, and the tp, ts, NM
    and cg tags."""
    out = []
    for line in text.splitlines():
        f = line.split("\t")
        tags = dict((t[:2], t[5:]) for t in f[12:])
        out.append(tuple(f[:12]) + (tags.get("tp"), tags.get("ts"), tags.get("NM"), tags.get("cg")))
    return sorted(out)


def check_paf(preset, ref, reads, ours, mate2=""):
    """'' when our PAF records equal minimap2 -c's for the same inputs."""
    cmd = [_minimap2(), "-c", "-x", preset, "-t", "1", ref, reads] + ([mate2] if mate2 else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return "minimap2 failed: " + r.stderr[-1000:]
    want, got = _paf(r.stdout), _paf(ours)
    if want == got:
        return ""
    w, g = set(want), set(got)
    return "%d records vs minimap2's %d; ours only %r; minimap2 only %r" % (len(got), len(want), sorted(g - w)[:2], sorted(w - g)[:2])


def spliced_count(paf):
    return sum(1 for l in paf.splitlines() if "N" in l.split("cg:Z:")[-1] and "cg:Z:" in l)


def fasta_lengths(path):
    names, lens = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
                lens.append(0)
            else:
                lens[-1] += len(line.strip())
    return (names, lens)
