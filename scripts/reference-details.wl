(* reference-details.wl -- curated detail layer for the API reference, merged
   on top of each symbol's ::usage by scripts/generate-reference.wls.

   Per symbol (all keys except "description" optional):
     "signatures" -> {call form, ...}        shown as a synopsis
     "description" -> markdown prose          (write it Markdown-safe)
     "args"        -> {{name, description}, ...}
     "options"     -> {{name, default, description}, ...}
     "messages"    -> {{"Symbol::tag", "when it is emitted"}, ...}
     "examples"    -> { "expr", {"expr", "display note"},
                        {"expr", "image", "path", "alt"}, ... }
                      a bare string is run and its output embedded; a
                      {expr, note} pair is run to validate but shows the note
                      instead of output (for front-end display forms); an
                      {expr, "image", path} triple embeds a committed image.
     "notes"       -> {note, ...}
     "seealso"     -> {SymbolName | {label, url}, ...}   (url for narrative pages)

   Every example is evaluated when the docs are built, so this file cannot
   drift from the code without the build failing. *)

<|

(* ============================= Algebras & generators ===================== *)

"LieAlgebra" -> <|
  "signatures" -> {"LieAlgebra[type, rank]"},
  "description" -> "The canonical constructor for a simple Lie algebra. `type` is one of the Cartan types and `rank` is a positive integer. Every algebra in the package normalizes to this form, so `SU[n]`, `SO[n]`, and `Sp[n]` are convenience wrappers that evaluate to a `LieAlgebra[...]`.",
  "args" -> {
    {"type", "the Cartan type: \"A\", \"B\", \"C\", or \"D\""},
    {"rank", "a positive integer, the rank of the algebra"}},
  "messages" -> {
    {"LieAlgebra::badtype", "type is not one of \"A\", \"B\", \"C\", \"D\""},
    {"LieAlgebra::badrank", "rank is not a valid rank for the given type"}},
  "examples" -> {"SU[3]", "LieAlgebra[\"B\", 2]"},
  "notes" -> {"\"A\" is su, \"B\" is so of odd dimension, \"C\" is sp, and \"D\" is so of even dimension."},
  "seealso" -> {"SU", "SO", "Sp", "Rank", "CartanMatrix", {"Concepts: the ABCD classification", "../concepts.md#simple-lie-algebras-and-the-abcd-classification"}}|>,

"SU" -> <|
  "signatures" -> {"SU[n]"},
  "description" -> "The special unitary algebra su(n), the type A(n-1) algebra. Shorthand that normalizes to `LieAlgebra[\"A\", n-1]`.",
  "args" -> {{"n", "an integer >= 2"}},
  "messages" -> {{"SU::baddim", "n is not an integer >= 2"}},
  "examples" -> {"SU[3]", "SU[3] === LieAlgebra[\"A\", 2]"},
  "seealso" -> {"LieAlgebra", "SO", "Sp", "Generators", "Irrep"}|>,

"SO" -> <|
  "signatures" -> {"SO[n]"},
  "description" -> "The special orthogonal algebra so(n). Covers both families: type B for odd n >= 3 and type D for even n >= 4.",
  "args" -> {{"n", "an integer >= 3"}},
  "messages" -> {{"SO::baddim", "n is not an integer >= 3"}},
  "examples" -> {"SO[5]", "SO[6]"},
  "notes" -> {"so(odd) and so(even) are different Cartan types (B and D), so their root systems and irreps differ."},
  "seealso" -> {"LieAlgebra", "SU", "Sp", "Generators", "Irrep"}|>,

"Sp" -> <|
  "signatures" -> {"Sp[n]"},
  "description" -> "The symplectic algebra sp(n), the type C(n/2) algebra. Defined for even n >= 2.",
  "args" -> {{"n", "an even integer >= 2 (the size of the defining representation)"}},
  "messages" -> {{"Sp::evenrank", "n is not an even integer >= 2"}},
  "examples" -> {"Sp[4]", "Sp[4] === LieAlgebra[\"C\", 2]"},
  "seealso" -> {"LieAlgebra", "SU", "SO", "Generators", "Irrep"}|>,

"Generators" -> <|
  "signatures" -> {"Generators[g]", "Generators[g, \"Standard\"]", "Generators[g, \"CartanWeyl\"]", "Generators[g, \"Chevalley\"]"},
  "description" -> "The generators of `g` in one of three bases. The scheme defaults to `\"Standard\"` (`Generators[g]` is `Generators[g, \"Standard\"]`): the defining-representation generators, as a list of exact matrices. `\"CartanWeyl\"` and `\"Chevalley\"` instead return an association with keys `\"Cartan\"`, `\"Raising\"`, `\"Lowering\"`.",
  "args" -> {
    {"g", "a Lie algebra, e.g. SU[3]"},
    {"scheme", "optional: \"Standard\" (default), \"CartanWeyl\", or \"Chevalley\""}},
  "options" -> {{"\"Realization\"", "\"Diagonal\"", "for so/sp in the **Chevalley** scheme, selects the Diagonal realization (diagonal Cartan generators) or the Antisymmetric one; ignored for the Standard/CartanWeyl schemes and for su"}},
  "messages" -> {
    {"Generators::badscheme", "the scheme is not \"Standard\", \"CartanWeyl\", or \"Chevalley\""},
    {"Generators::badrealization", "\"Realization\" is not \"Diagonal\" or \"Antisymmetric\" (Chevalley scheme of so/sp)"}},
  "examples" -> {"Generators[SU[2]]", "Keys[Generators[SU[3], \"Chevalley\"]]", "Generators[SU[3], \"Standard\"] === Generators[SU[3]]"},
  "notes" -> {
    "The \"Realization\" option must be given with an explicit scheme, e.g. Generators[g, \"Chevalley\", \"Realization\" -> \"Antisymmetric\"]; passing it with no scheme is read as the scheme argument (a badscheme error). BasisTransform[g] relates the two Chevalley realizations.",
    "The CartanWeyl and Chevalley shorthands take no options; use Generators[g, ...] directly to pass \"Realization\"."},
  "seealso" -> {"CartanWeyl", "Chevalley", "BasisTransform", "RepresentationMatrices", {"Concepts: the three bases", "../concepts.md#the-three-generator-bases"}, {"Tutorial: getting started", "../tutorials/getting-started.md"}}|>,

"CartanWeyl" -> <|
  "signatures" -> {"CartanWeyl[g]"},
  "description" -> "Shorthand for `Generators[g, \"CartanWeyl\"]`. Returns the Cartan-Weyl basis as an association: the diagonal Cartan generators (`\"Cartan\"`), the raising operators (`\"Raising\"`), and the lowering operators (`\"Lowering\"`).",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"Keys[CartanWeyl[SU[2]]]", "CartanWeyl[SU[2]][\"Cartan\"]"},
  "notes" -> {"A single-argument shorthand: it takes no options. To choose a realization, call Generators[g, \"CartanWeyl\", ...] directly."},
  "seealso" -> {"Generators", "Chevalley", {"Guided tour: Cartan-Weyl basis", "../walkthrough.md#cartanweyl-basis"}}|>,

"Chevalley" -> <|
  "signatures" -> {"Chevalley[g]"},
  "description" -> "Shorthand for `Generators[g, \"Chevalley\"]`. Returns the Chevalley basis as an association of `\"Cartan\"`, `\"Raising\"`, `\"Lowering\"` generators, one triple per simple root. This is the basis the representation engine uses; the default (Diagonal) realization has diagonal Cartan generators.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"Chevalley[SU[2]]", "DiagonalMatrixQ /@ Chevalley[SO[5]][\"Cartan\"]"},
  "notes" -> {"A single-argument shorthand: it takes no options. For the Antisymmetric realization of so/sp, call Generators[g, \"Chevalley\", \"Realization\" -> \"Antisymmetric\"]; BasisTransform[g] maps it to the Diagonal one."},
  "seealso" -> {"Generators", "CartanWeyl", "RepresentationMatrices", {"Guided tour: Chevalley basis", "../walkthrough.md#chevalley-basis"}}|>,

"BasisTransform" -> <|
  "signatures" -> {"BasisTransform[g]"},
  "description" -> "The change-of-basis matrix conjugating the Antisymmetric realization of an so/sp Chevalley basis into the Diagonal one (`u.a.Inverse[u]` carries one to the other). For su it is the identity.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"BasisTransform[SU[2]]", "Dimensions[BasisTransform[SO[5]]]"},
  "seealso" -> {"Generators", "Chevalley"}|>,

(* ================================ Root system ============================ *)

"Rank" -> <|
  "signatures" -> {"Rank[g]"},
  "description" -> "The rank of `g` -- the dimension of its Cartan subalgebra, equal to the number of simple roots.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"Rank[SO[5]]", "Rank[SU[4]]"},
  "seealso" -> {"LieAlgebraDimension", "CartanMatrix", "SimpleRoots"}|>,

"LieAlgebraDimension" -> <|
  "signatures" -> {"LieAlgebraDimension[g]"},
  "description" -> "The dimension of `g` as a vector space -- equivalently, the number of generators.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"LieAlgebraDimension[SU[3]]", "LieAlgebraDimension[SO[5]]"},
  "notes" -> {"Equals the number of defining-representation generators returned by Generators[g], and (# positive roots) * 2 + rank."},
  "seealso" -> {"Rank", "Generators", "PositiveRoots"}|>,

"CartanMatrix" -> <|
  "signatures" -> {"CartanMatrix[g]"},
  "description" -> "The Cartan matrix of `g`: the rank-by-rank integer matrix `A[[i,j]] = 2 (a_i, a_j) / (a_j, a_j)` of pairings between simple roots and coroots.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"CartanMatrix[SU[3]]", "CartanMatrix[SO[5]]"},
  "notes" -> {"For the non-simply-laced types (B, C) this is the transpose of the matrix some references (e.g. LieART) print -- CartanMatrix[SO[5]] is {{2, -2}, {-1, 2}} -- equivalently the package uses [H_i, E_j] = A[[j,i]] E_j. It is internally consistent in this convention."},
  "seealso" -> {"SimpleRoots", "Rank", {"Concepts: conventions", "../concepts.md#conventions"}}|>,

"SimpleRoots" -> <|
  "signatures" -> {"SimpleRoots[g]"},
  "description" -> "The simple roots of `g`, as vectors in the Euclidean (orthonormal) basis. Ordered by the Bourbaki labelling of the Dynkin diagram.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"SimpleRoots[SU[3]]", "SimpleRoots[Sp[4]]"},
  "seealso" -> {"PositiveRoots", "CartanMatrix", "FundamentalWeights"}|>,

"PositiveRoots" -> <|
  "signatures" -> {"PositiveRoots[g]"},
  "description" -> "All positive roots of `g`, as vectors in the Euclidean basis. There are (dim(g) - rank(g))/2 of them.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"PositiveRoots[SU[3]]", "PositiveRoots[SO[5]]"},
  "seealso" -> {"SimpleRoots", "FundamentalWeights"}|>,

"FundamentalWeights" -> <|
  "signatures" -> {"FundamentalWeights[g]"},
  "description" -> "The fundamental weights of `g`, as vectors in the Euclidean basis. They are dual to the simple coroots, and the Dynkin labels of a weight are its coordinates in this basis.",
  "args" -> {{"g", "a Lie algebra"}},
  "examples" -> {"FundamentalWeights[SU[3]]", "FundamentalWeights[SO[5]]"},
  "notes" -> {"Half-integer entries (e.g. {1/2, 1/2} for SO[5]) are the spinor fundamental weights."},
  "seealso" -> {"SimpleRoots", "Irrep"}|>,

(* ============================== Representations =========================== *)

"Irrep" -> <|
  "signatures" -> {"Irrep[g, w]"},
  "description" -> "Represents the irreducible representation of `g` with highest weight `w`. This is an inert label: a valid `Irrep[...]` returns unchanged, and the data functions (`RepresentationDimension`, `WeightSystem`, `CasimirEigenvalue`, `RepresentationMatrices`) take it and compute from it.",
  "args" -> {
    {"g", "a Lie algebra"},
    {"w", "the highest weight as Dynkin labels: a list of Rank[g] non-negative integers"}},
  "messages" -> {{"Irrep::badweight", "w is not a list of Rank[g] non-negative integers (wrong length or a negative entry)"}},
  "examples" -> {"Irrep[SU[3], {1, 1}]", "RepresentationDimension[Irrep[SO[6], {1, 0, 0}]]"},
  "notes" -> {"The highest weight must be non-negative; the weights returned by WeightSystem may have negative entries."},
  "seealso" -> {"RepresentationDimension", "WeightSystem", "CasimirEigenvalue", "RepresentationMatrices", "HighestWeight", {"Tutorial: representations", "../tutorials/representations.md"}}|>,

"HighestWeight" -> <|
  "signatures" -> {"HighestWeight[Irrep[g, w]]"},
  "description" -> "Returns the highest weight `w` (Dynkin labels) of an irrep.",
  "args" -> {{"ir", "an Irrep[g, w]"}},
  "examples" -> {"HighestWeight[Irrep[SU[3], {1, 1}]]"},
  "seealso" -> {"Irrep"}|>,

"RepresentationDimension" -> <|
  "signatures" -> {"RepresentationDimension[Irrep[g, w]]"},
  "description" -> "The dimension of the irrep, computed in closed form from the Weyl dimension formula (no enumeration of states needed).",
  "args" -> {{"ir", "an Irrep[g, w]"}},
  "examples" -> {"RepresentationDimension[Irrep[SU[3], {1, 1}]]", "RepresentationDimension[Irrep[SO[5], {0, 1}]]"},
  "notes" -> {"Reaches the spinor representations: the so(5) irrep {0, 1} above is the 4-dimensional spinor."},
  "seealso" -> {"WeightSystem", "CasimirEigenvalue", "Irrep", {"Tutorial: representations", "../tutorials/representations.md"}}|>,

"RepresentationMatrices" -> <|
  "signatures" -> {"RepresentationMatrices[Irrep[g, w]]"},
  "description" -> "The Chevalley generators of `g` (one `H`, `E`, `F` per simple root) as explicit matrices in the irrep, returned as an association with keys `\"Cartan\"`, `\"Raising\"`, `\"Lowering\"`. The default basis is orthonormal, with each `E` the conjugate transpose of the corresponding `F` and each `H` diagonal; each value is a list of `Rank[g]` matrices of size `RepresentationDimension` square.",
  "args" -> {{"ir", "an Irrep[g, w]"}},
  "examples" -> {"RepresentationMatrices[Irrep[SU[2], {1}]]", "Diagonal[RepresentationMatrices[Irrep[SO[5], {0, 1}]][\"Cartan\"][[1]]]"},
  "notes" -> {
    "Built by the abstract highest-weight (Shapovalov) construction, so it works for every classical type and every irrep, spinors included -- the so(5) spinor {0, 1} above is a genuine 4x4 example.",
    "Cost grows quickly with the dimension of the irrep; build small irreps, or wrap larger ones in TimeConstrained."},
  "seealso" -> {"Chevalley", "WeightSystem", "CasimirEigenvalue", "Irrep", {"Guided tour", "../walkthrough.md#representations"}, {"Tutorial: physics applications", "../tutorials/physics-applications.md"}}|>,

"CasimirEigenvalue" -> <|
  "signatures" -> {"CasimirEigenvalue[Irrep[g, w]]"},
  "description" -> "The eigenvalue of the quadratic Casimir operator on the irrep, normalized so that long roots have squared length 2 (the mathematicians' normalization).",
  "args" -> {{"ir", "an Irrep[g, w]"}},
  "examples" -> {"CasimirEigenvalue[Irrep[SU[3], {1, 1}]]", "CasimirEigenvalue[Irrep[SU[2], {1}]]", "CasimirEigenvalue[Irrep[SO[5], {0, 1}]]"},
  "notes" -> {"With this normalization the su(2) spin-1/2 value is 3/2 -- twice the physics value 3/4."},
  "seealso" -> {"RepresentationDimension", "WeightSystem", "Irrep"}|>,

"WeightSystem" -> <|
  "signatures" -> {"WeightSystem[Irrep[g, w]]"},
  "description" -> "The complete weight system of the irrep: an association from each weight (a length-`Rank[g]` Dynkin-label list) to its multiplicity (a positive integer), obtained from the Freudenthal multiplicity recursion. The multiplicities sum to `RepresentationDimension`.",
  "args" -> {{"ir", "an Irrep[g, w]"}},
  "examples" -> {"WeightSystem[Irrep[SU[3], {1, 0}]]", "WeightSystem[Irrep[SU[3], {1, 1}]]", "WeightSystem[Irrep[SO[5], {0, 1}]]"},
  "notes" -> {"A weight may have multiplicity above 1: in the adjoint {1, 1} of su(3) the zero weight {0, 0} appears twice."},
  "seealso" -> {"RepresentationDimension", "HighestWeight", "Irrep", {"Tutorial: physics applications", "../tutorials/physics-applications.md"}}|>,

(* =============================== Young tableaux ========================== *)

"Tableau" -> <|
  "signatures" -> {"Tableau[rows]"},
  "description" -> "A Young tableau given as a list of row-lists. Used to specify the Young symmetrizer applied by `TableauPermute`.",
  "args" -> {{"rows", "a list of rows, e.g. {{1, 2}, {3}}"}},
  "examples" -> {"Tableau[{{1, 2}, {3}}]", "TableauPermute[Tableau[{{1, 2}, {3}}], Psi[1, 2, 3]]"},
  "seealso" -> {"TableauPermute", "TensorTableau", {"Guided tour", "../walkthrough.md#young-tableaux-and-wave-functions"}}|>,

"TensorTableau" -> <|
  "signatures" -> {"TensorTableau[rows]"},
  "description" -> "A tensor-product Young tableau with explicit index entries -- a compact notation for a many-body wave function that also displays its permutation symmetry. Linear combinations of `TensorTableau` expressions are the working objects of the toolkit.",
  "args" -> {{"rows", "a list of rows of index entries, e.g. {{1, 2}, {3}}"}},
  "examples" -> {"TensorTableau[{{1, 2}, {3}}]", "ToTensor[TensorTableau[{{1, 2}, {3}}]]"},
  "seealso" -> {"ToTensor", "TableauForm", "TableauDot", "TableauNormalization", {"Tutorial: Young tableaux", "../tutorials/young-tableaux.md"}}|>,

"TableauForm" -> <|
  "signatures" -> {"TableauForm[t]"},
  "description" -> "Displays a `TensorTableau` (or a linear combination of them) as a labelled Young-tableau grid. This is a display form; use `ToTensor` to get the same state as data.",
  "args" -> {{"t", "a TensorTableau or linear combination"}},
  "examples" -> {
    {"TableauForm[TensorTableau[{{1, 2}, {3}}]]", "renders as a labelled Young-tableau grid in the notebook front end (top row 1, 2; second row 3). The same state, as data:"},
    "ToTensor[TensorTableau[{{1, 2}, {3}}]]"},
  "seealso" -> {"TensorTableau", "ToTensor"}|>,

"TableauPermute" -> <|
  "signatures" -> {"TableauPermute[t, v]"},
  "description" -> "Applies the Young symmetrizer of the tableau `t` to a basis tensor `v`, giving the (unnormalized) wave function as a combination of `Psi` tensors. A repeated entry within a column makes the result vanish, since columns are antisymmetrized.",
  "args" -> {
    {"t", "a Tableau giving the symmetry"},
    {"v", "a Psi basis tensor"}},
  "examples" -> {"TableauPermute[Tableau[{{1, 2}, {3}}], Psi[1, 1, 2]]"},
  "seealso" -> {"Tableau", "Psi", "ToTensor"}|>,

"TableauNormalization" -> <|
  "signatures" -> {"TableauNormalization[t]"},
  "description" -> "Normalizes a `TensorTableau` expression to unit norm.",
  "args" -> {{"t", "a TensorTableau or linear combination"}},
  "examples" -> {"TableauNormalization[TensorTableau[{{1, 2}, {2}}]]"},
  "seealso" -> {"TensorNorm", "TableauDot", "TableauOrthogonalization"}|>,

"TableauOrthogonalization" -> <|
  "signatures" -> {"TableauOrthogonalization[t1, t2]"},
  "description" -> "One Gram-Schmidt step: returns `{t1, t2 - <t2,t1>/<t1,t1> t1}`, keeping `t1` and orthogonalizing `t2` against it. Used to pick a basis for a degenerate (repeated) weight space.",
  "args" -> {
    {"t1", "a TensorTableau expression to keep"},
    {"t2", "a TensorTableau expression to orthogonalize against t1"}},
  "examples" -> {"TableauOrthogonalization[2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}], TensorTableau[{{1, 3}, {2}}] + TensorTableau[{{1, 2}, {3}}]]"},
  "seealso" -> {"TableauDot", "TableauNormalization", {"Guided tour: orthogonalizing a degenerate weight", "../walkthrough.md#the-11-representation-of-su3"}}|>,

"TableauDot" -> <|
  "signatures" -> {"TableauDot[t1, t2]"},
  "description" -> "The inner product of two `TensorTableau` expressions.",
  "args" -> {
    {"t1", "a TensorTableau expression"},
    {"t2", "a TensorTableau expression"}},
  "examples" -> {"TableauDot[TensorTableau[{{1, 2}, {3}}], TensorTableau[{{1, 2}, {3}}]]"},
  "seealso" -> {"TensorDot", "TableauNormalization", "TableauOrthogonalization"}|>,

"ToTensor" -> <|
  "signatures" -> {"ToTensor[t]"},
  "description" -> "Converts a `TensorTableau` (or linear combination) into an explicit linear combination of `Psi` basis tensors.",
  "args" -> {{"t", "a TensorTableau or linear combination"}},
  "examples" -> {"ToTensor[TensorTableau[{{1, 2}, {3}}]]"},
  "seealso" -> {"TensorTableau", "Psi", "TableauPermute"}|>,

"TensorDot" -> <|
  "signatures" -> {"TensorDot[p1, p2]"},
  "description" -> "The Hermitian inner product of two linear combinations of `Psi` tensors. Distinct basis tensors are orthonormal.",
  "args" -> {
    {"p1", "a linear combination of Psi tensors"},
    {"p2", "a linear combination of Psi tensors"}},
  "examples" -> {"TensorDot[Psi[1, 2, 3], Psi[1, 2, 3]]", "TensorDot[Psi[1, 2], Psi[2, 1]]"},
  "seealso" -> {"TableauDot", "TensorNorm", "Psi"}|>,

"TensorNorm" -> <|
  "signatures" -> {"TensorNorm[p]"},
  "description" -> "The norm of a linear combination of `Psi` tensors, i.e. Sqrt[TensorDot[p, p]].",
  "args" -> {{"p", "a linear combination of Psi tensors"}},
  "examples" -> {"TensorNorm[Psi[1, 2] + Psi[2, 1]]"},
  "seealso" -> {"TensorDot", "TableauNormalization"}|>,

"Psi" -> <|
  "signatures" -> {"Psi[i1, i2, ...]"},
  "description" -> "A basis tensor |i1, i2, ...> in the tensor-product space -- the elementary building block of the many-body wave functions. Distinct `Psi` tensors are orthonormal under `TensorDot`.",
  "args" -> {{"i1, i2, ...", "the indices labelling the basis tensor"}},
  "examples" -> {"Psi[1, 2, 3]", "TensorDot[Psi[1, 2], Psi[1, 2]]"},
  "seealso" -> {"ToTensor", "TensorDot", "TensorNorm", "TableauPermute"}|>

|>
