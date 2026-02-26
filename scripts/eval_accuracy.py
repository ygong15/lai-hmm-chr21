import argparse
import csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="data/sim_truth.tsv")
    ap.add_argument("--pred", default="data/pred.tsv")
    args = ap.parse_args()

    # Read truth into dict keyed by (CHR, POS)
    truth = {}
    with open(args.truth, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = (row["CHR"], row["POS"])
            truth[key] = (row["TRUE1"], row["TRUE2"])

    n = 0
    c1 = 0
    c2 = 0

    # confusion counts per hap: (true,pred)->count
    conf1 = {}
    conf2 = {}

    with open(args.pred, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = (row["CHR"], row["POS"])
            if key not in truth:
                continue
            t1, t2 = truth[key]
            p1, p2 = row["POP1"], row["POP2"]

            n += 1
            c1 += (p1 == t1)
            c2 += (p2 == t2)

            conf1[(t1, p1)] = conf1.get((t1, p1), 0) + 1
            conf2[(t2, p2)] = conf2.get((t2, p2), 0) + 1

    if n == 0:
        raise SystemExit("No overlapping sites between truth and pred.")

    acc1 = c1 / n
    acc2 = c2 / n
    mean = (acc1 + acc2) / 2

    print(f"Sites compared: {n}")
    print(f"Accuracy POP1: {acc1:.6f}")
    print(f"Accuracy POP2: {acc2:.6f}")
    print(f"Mean accuracy: {mean:.6f}")

    def show_conf(conf, label):
        # order rows for readability
        pairs = [("A","A"),("A","B"),("B","A"),("B","B")]
        print(f"\nConfusion {label} (TRUE,PRED):")
        for t,p in pairs:
            print(f"  ({t},{p}) {conf.get((t,p),0)}")

    show_conf(conf1, "POP1")
    show_conf(conf2, "POP2")

if __name__ == "__main__":
    main()