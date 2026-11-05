"""Command line: ``python -m morpho_gruson.cli run --raw data/raw --config config/pipeline.yaml --out data/processed``"""
import argparse, warnings
from .pipeline import run
from .export import export_all


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["run"])
    ap.add_argument("--raw", default="data/raw"); ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--out", default="data/processed")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    R = run(a.config, a.raw)
    nb, nh = export_all(R, a.out)
    print(f"{len(R['summary'])} species | {nb} brightness sweeps | {nh} hue sweeps -> {a.out}")


if __name__ == "__main__":
    main()
