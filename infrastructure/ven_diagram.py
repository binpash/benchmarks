import csv
import matplotlib.pyplot as plt
import sys
import networkx as nx
import math

benchmark_csv = sys.argv[1] if len(sys.argv) > 1 else None
out_filename = "category_venn_diagram.pdf"
assert benchmark_csv, "Please provide the path to the benchmarks CSV file as an argument."

# 1) Define categories and accumulate total LoC + co-occurrences
categories = ["DA", "AN", "ML", "TP", "SA", "CI", "MI"]
loc_sums = {cat: 0 for cat in categories}
cooccurrences = set()

with open(benchmark_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        first, second = row['first_class'], row['secondary_class']
        loc = int(row['loc'])
        if first in loc_sums:   loc_sums[first]  += loc
        if second in loc_sums:  loc_sums[second] += loc
        if first in categories and second in categories and first != second:
            cooccurrences.add(tuple(sorted([first, second])))

# 2) Build graph of categories
G = nx.Graph()
for cat in categories:
    G.add_node(cat, loc=loc_sums[cat])

for a, b in cooccurrences:
    G.add_edge(a, b)

# 3) Initial layout
pos = nx.spring_layout(G, seed=42)

# 4) Compute circle radii (area ∝ LoC)
max_loc = max(loc_sums.values())
max_display_radius = 0.15          # largest circle ≈15% of figure width
scale = math.sqrt(max_loc / math.pi) / max_display_radius
radii = {n: math.sqrt(data['loc']/math.pi)/scale for n, data in G.nodes(data=True)}

# 5) Nudge overlapping nodes together
#    For each edge, if dist > r1 + r2, move them closer
for _ in range(100):
    for u, v in G.edges():
        x1, y1 = pos[u];  x2, y2 = pos[v]
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        desired = (radii[u] + radii[v]) * 0.9  # 90% overlap factor
        if dist > desired and dist > 1e-6:
            # shift each node halfway toward one another
            shift = 0.5 * (dist - desired)
            ux, uy = dx/dist, dy/dist
            pos[u] = (x1 + ux*shift, y1 + uy*shift)
            pos[v] = (x2 - ux*shift, y2 - uy*shift)

# 6) Draw
fig, ax = plt.subplots(figsize=(10,10))
all_x, all_y = [], []

for node in G.nodes():
    x, y = pos[node]
    r = radii[node]
    circle = plt.Circle((x, y), r, alpha=0.5, edgecolor='black')
    ax.add_patch(circle)
    ax.text(x, y, node, ha='center', va='center', fontsize=12, weight='bold')
    all_x += [x-r, x+r]
    all_y += [y-r, y+r]

# expand limits so no circle is cut off
margin = 0.05 * max(max(all_x)-min(all_x), max(all_y)-min(all_y))
ax.set_xlim(min(all_x)-margin, max(all_x)+margin)
ax.set_ylim(min(all_y)-margin, max(all_y)+margin)

ax.set_aspect('equal')
ax.axis('off')
plt.tight_layout()
plt.savefig(out_filename)
print(f"Wrote category‐based overlapping diagram to {out_filename}")
