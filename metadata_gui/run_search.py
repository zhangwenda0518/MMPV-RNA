"""Run gsa_sra.search.py and capture output."""
import subprocess, sys, os, time

script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "public_metadata_pipeline", "gsa_sra.search.py")
outdir = os.path.join(os.path.dirname(script), "test_output")
os.makedirs(outdir, exist_ok=True)
log_path = os.path.join(outdir, "run.log")

with open(log_path, "w", encoding="utf-8") as f:
    f.write(f"=== START {time.ctime()} ===\n")
    f.write(f"Script: {script}\n")
    f.write(f"Outdir: {outdir}\n")
    f.flush()

    p = subprocess.Popen(
        [sys.executable, "-u", script, "-q", "Lycium chinense", "--db", "both", "-o", outdir],
        stdout=f, stderr=subprocess.STDOUT,
        cwd=os.path.dirname(script))

    f.write(f"PID: {p.pid}\n")
    f.flush()

    try:
        rc = p.wait(timeout=90)
        f.write(f"\n=== DONE {time.ctime()} RC={rc} ===\n")
    except subprocess.TimeoutExpired:
        p.kill()
        f.write(f"\n=== TIMEOUT {time.ctime()} ===\n")

    # List output files
    f.write("\n--- Files ---\n")
    for fn in os.listdir(outdir):
        fp = os.path.join(outdir, fn)
        f.write(f"  {fn} ({os.path.getsize(fp)} bytes)\n")

print("Done. Log:", log_path)
