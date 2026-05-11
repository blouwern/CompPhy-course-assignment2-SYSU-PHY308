import argparse
import csv
import math
import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def load_observable(path: str) -> np.ndarray:
    values = []
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            values.append(float(row["observable"]))
    return np.asarray(values, dtype=float)


def autocorrelation(x: np.ndarray, max_lag: int) -> np.ndarray:
    x = x - np.mean(x)
    n = len(x)
    if n == 0:
        return np.zeros(1)
    var = np.var(x)
    if var == 0:
        return np.zeros(max_lag + 1)
    acf = np.zeros(max_lag + 1)
    acf[0] = 1.0
    for lag in range(1, max_lag + 1):
        acf[lag] = np.dot(x[:-lag], x[lag:]) / ((n - lag) * var)
    return acf


def integrated_autocorr_time(acf: np.ndarray) -> float:
    if len(acf) <= 1:
        return 1.0
    tau = 1.0
    for lag in range(1, len(acf)):
        if acf[lag] <= 0.0:
            break
        tau += 2.0 * acf[lag]
    return tau


def analyze_file(path: str, max_lag: int, out_dir: str) -> Tuple[float, int]:
    series = load_observable(path)
    if len(series) == 0:
        return 1.0, 1
    max_lag = min(max_lag, len(series) - 1)
    acf = autocorrelation(series, max_lag)
    tau = integrated_autocorr_time(acf)
    thin = max(1, int(math.ceil(2.0 * tau)))

    base = os.path.splitext(os.path.basename(path))[0]

    plt.figure(figsize=(8, 3))
    plt.plot(series[: min(2000, len(series))], lw=0.8)
    plt.title(f"Trace: {base}")
    plt.xlabel("Sample")
    plt.ylabel("Observable")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"trace_{base}.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(6, 3))
    plt.plot(acf, lw=1.0)
    plt.title(f"ACF: {base}")
    plt.xlabel("Lag")
    plt.ylabel("ACF")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"acf_{base}.png"), dpi=150)
    plt.close()

    return tau, thin


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", help="CSV files from the sampler")
    parser.add_argument("--max-lag", type=int, default=2000)
    parser.add_argument("--out-dir", default="figures")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    rows: List[Tuple[str, float, int]] = []
    for path in args.inputs:
        tau, thin = analyze_file(path, args.max_lag, args.out_dir)
        rows.append((path, tau, thin))

    print("file,tau_int,thin")
    for path, tau, thin in rows:
        print(f"{path},{tau:.3f},{thin}")


if __name__ == "__main__":
    main()
