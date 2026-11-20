# The Complexity of Minimum-Energy Strong Orientations

Yitzchak Shmalo

For an orientation of a graph, let $f$ solve $\Delta f=d^{\mathrm{out}}-d^{\mathrm{in}}$; the energy of the orientation is
$\langle f,\Delta f\rangle$. Amini showed that the squared covering radius of the lattice of integer flows is an affine
function of the least energy of a strongly connected orientation, and asked whether that least energy can be computed in
polynomial time. Deciding whether it is at most a given rational number is NP-complete, even on simple bipartite planar
graphs that become a tree after one vertex is deleted, so deciding whether the squared covering radius of the lattice of
integer flows (or of integer cuts) is at most a given number is coNP-complete. At bounded treewidth the least energy has a
fully polynomial-time approximation scheme.

`main.tex` and `main.pdf` are the paper. `check_triangle_ring.py` uses only Python's standard library and computes the
values of Example 3.4 from the definition of the energy, in exact rational arithmetic:

    python check_triangle_ring.py

Its output is recorded in `check_output.txt`.
