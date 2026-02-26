import argparse
import csv
import matplotlib.pyplot as plt

COLOR = {"A": "tab:blue", "B": "tab:red"}

def read_truth(truth_fp):
    rows = []
    with open(truth_fp, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pos = int(row["POS"])
            rows.append((pos, row["TRUE1"], row["TRUE2"]))
    return rows

def read_pred(pred_fp):
    rows = []
    with open(pred_fp, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pos = int(row["POS"])
            rows.append((pos, row["POP1"], row["POP2"]))
    return rows

def compress_segments(pos_state_list):
    if not pos_state_list:
        return []
    segs = []
    cur_state = pos_state_list[0][1]
    start = pos_state_list[0][0]
    prev = pos_state_list[0][0]
    for pos, state in pos_state_list[1:]:
        if state != cur_state:
            segs.append((start, prev, cur_state))
            cur_state = state
            start = pos
        prev = pos
    segs.append((start, prev, cur_state))
    return segs

def draw_segments(ax, segs, y, mb):
    for start, end, state in segs:
        x0, x1 = start, end
        if mb:
            x0 /= 1e6
            x1 /= 1e6
        ax.plot([x0, x1], [y, y],
                linewidth=10,
                solid_capstyle="butt",
                color=COLOR.get(state, "black"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="data/sim_truth.tsv")
    ap.add_argument("--pred", default="data/sim_pred.tsv")
    ap.add_argument("--out", default="results/truth_vs_pred.png")
    ap.add_argument("--title", default="Truth vs Pred (A=AFR, B=EUR)")
    ap.add_argument("--mb", action="store_true")
    args = ap.parse_args()

    truth = read_truth(args.truth)
    pred = read_pred(args.pred)

    if len(truth) != len(pred):
        raise SystemExit(f"Length mismatch truth={len(truth)} pred={len(pred)}")

    for (pt, _, _), (pp, _, _) in zip(truth, pred):
        if pt != pp:
            raise SystemExit(f"POS mismatch truth={pt} pred={pp}")

    pos = [p for (p, _, _) in truth]

    true1 = [(pos[i], truth[i][1]) for i in range(len(pos))]
    pred1 = [(pos[i], pred[i][1]) for i in range(len(pos))]
    true2 = [(pos[i], truth[i][2]) for i in range(len(pos))]
    pred2 = [(pos[i], pred[i][2]) for i in range(len(pos))]

    seg_true1 = compress_segments(true1)
    seg_pred1 = compress_segments(pred1)
    seg_true2 = compress_segments(true2)
    seg_pred2 = compress_segments(pred2)

    fig, ax = plt.subplots()
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["True Hap1", "Pred Hap1", "True Hap2", "Pred Hap2"])

    draw_segments(ax, seg_true1, 0, args.mb)
    draw_segments(ax, seg_pred1, 1, args.mb)
    draw_segments(ax, seg_true2, 2, args.mb)
    draw_segments(ax, seg_pred2, 3, args.mb)

    xmin = pos[0] / (1e6 if args.mb else 1.0)
    xmax = pos[-1] / (1e6 if args.mb else 1.0)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Position (Mb)" if args.mb else "Position (bp)")
    ax.set_title(args.title)

    # Manual legend (matches COLOR map)
    ax.plot([], [], color=COLOR["A"], linewidth=10, label="A")
    ax.plot([], [], color=COLOR["B"], linewidth=10, label="B")
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Wrote plot: {args.out}")

if __name__ == "__main__":
    main()