 
import argparse
import glob
from collections import Counter
from ete3 import Tree
 
 
def canonical_subtree(node):
    """Return a sorted topology string for any rooted subtree."""
    if node.is_leaf():
        return node.name
 
    children = [canonical_subtree(ch) for ch in node.children]
    children = sorted(children)
    return "(" + ",".join(children) + ")"
 
 
def classify(tree, outgroup):
    tree.set_outgroup(outgroup)
 
    leaves = sorted(tree.get_leaf_names())
 
    if outgroup not in leaves:
        raise ValueError(f"Outgroup {outgroup} not found in tree: {leaves}")
 
    ingroup = [t for t in leaves if t != outgroup]
 
    if len(ingroup) != 4:
        raise ValueError(f"Tree does not contain exactly 4 ingroup taxa: {leaves}")
 
    # After rooting on outgroup, the non-outgroup child of the root is the ingroup topology
    root_children = tree.children
    ingroup_nodes = [
        ch for ch in root_children
        if outgroup not in ch.get_leaf_names()
    ]
 
    if len(ingroup_nodes) != 1:
        raise ValueError(f"Could not identify single ingroup clade after rooting: {leaves}")
 
    ingroup_topo = canonical_subtree(ingroup_nodes[0])
 
    return f"{outgroup},{ingroup_topo}"
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Classify rooted 5-taxon gene trees and output topology counts."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing .treefile gene trees"
    )
    parser.add_argument(
        "-o", "--outgroup",
        required=True,
        help="Name of the outgroup taxon, e.g. Ficus"
    )
    parser.add_argument(
        "--output",
        default="topology_counts.tsv",
        help="Output TSV path"
    )
 
    args = parser.parse_args()
 
    pattern = f"{args.input.rstrip('/')}/*.treefile"
    files = glob.glob(pattern)
 
    if not files:
        print(f"No .treefile files found in {args.input}")
        with open(args.output, "w") as outfh:
            outfh.write("topology\tcount\n")
        return
 
    counts = Counter()
 
    for f in files:
        try:
            t = Tree(f)
            topo = classify(t, args.outgroup)
            counts[topo] += 1
        except Exception as e:
            print(f"Warning: failed to parse {f}: {e}")
 
    print("\nTopology counts:")
    for topo in sorted(counts.keys()):
        print(f"{topo}: {counts[topo]}")
 
    total = sum(counts.values())
    print(f"\nTotal trees processed: {total}")
 
    with open(args.output, "w") as outfh:
        outfh.write("topology\tcount\n")
        for topo in sorted(counts.keys()):
            outfh.write(f"{topo}\t{counts[topo]}\n")
 
    print(f"\nWrote TSV: {args.output}")
 
 
if __name__ == "__main__":
    main()
