# The Complexity of Minimum-Energy Strong Orientations

Yitzchak Shmalo

For an orientation of a finite graph, let $\Delta$ be the graph Laplacian and let $f$ solve $\Delta
f=d^{\mathrm{out}}-d^{\mathrm{in}}$, the out-degree minus the in-degree, and call $\langle f,\Delta f\rangle$ the energy
of the orientation. Amini showed that the covering radius of the lattice of integer flows is determined by the least
energy of a strongly connected orientation, and wrote that the complexity of finding such an orientation is unknown and
that he suspected it to be NP-hard. We show that the decision version is NP-complete, even for simple bipartite planar
graphs that become a tree after one vertex is deleted. Hence deciding whether the squared covering radius of the lattice
of integer flows is at most a given rational number is coNP-complete, and the same holds for lattices of integer cuts.
On graphs of bounded treewidth the least energy has a fully polynomial-time approximation scheme.

`main.tex` and `main.pdf` are the paper. `check_triangle_ring.py` uses only Python's standard library and computes the
values of Example 3.4 from the definition of the energy, in exact rational arithmetic:

    python check_triangle_ring.py

Its output is recorded in `check_output.txt`.
