import csv
import argparse


def read_truth(path):
    pos, h1, h2 = [], [], []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pos.append(int(row["POS"]))
            h1.append(row["TRUE1"])
            h2.append(row["TRUE2"])
    return pos, h1, h2


def read_pred(path):
    pos, h1, h2 = [], [], []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pos.append(int(row["POS"]))
            h1.append(row["POP1"])
            h2.append(row["POP2"])
    return pos, h1, h2


def count_switches(states):
    switches = 0
    for i in range(1, len(states)):
        if states[i] != states[i - 1]:
            switches += 1
    return switches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--pred", required=True)
    args = ap.parse_args()

    pos_t, t1, t2 = read_truth(args.truth)
    pos_p, p1, p2 = read_pred(args.pred)

    if pos_t != pos_p:
        raise ValueError("Positions do not match between truth and prediction.")

    n = len(pos_t)

    # Accuracy
    correct = 0
    for i in range(n):
        if t1[i] == p1[i]:
            correct += 1
        if t2[i] == p2[i]:
            correct += 1
    accuracy = correct / (2 * n)

    # Switch counts
    true_switches = count_switches(t1) + count_switches(t2)
    pred_switches = count_switches(p1) + count_switches(p2)

    skew = pred_switches / true_switches if true_switches > 0 else float("inf")

    print(f"Accuracy: {accuracy:.6f}")
    print(f"True switches: {true_switches}")
    print(f"Pred switches: {pred_switches}")
    print(f"Skew (pred/true): {skew:.6f}")


if __name__ == "__main__":
    main()