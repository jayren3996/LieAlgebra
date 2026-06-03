Needs["ClassicalLieAlgebra`"];

VerificationTest[ TableauPermute[Tableau[{{1, 2}, {3}}], Psi[1, 1, 2]],
  2 Psi[1, 1, 2] - Psi[2, 1, 1] - Psi[1, 2, 1], TestID -> "permute-112" ];
VerificationTest[ TensorNorm[Psi[1, 2, 3] - 2 Psi[3, 2, 1] + Psi[4, 5, 6] - 2 Psi[3, 2, 1]],
  Sqrt[18], TestID -> "tensornorm" ];
VerificationTest[ ToTensor[TensorTableau[{{1, 1}, {2}}]],
  2 Psi[1, 1, 2] - Psi[1, 2, 1] - Psi[2, 1, 1], TestID -> "totensor-11hw" ];
VerificationTest[
  TableauOrthogonalization[
    2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}],
    TensorTableau[{{1, 3}, {2}}] + TensorTableau[{{1, 2}, {3}}]],
  {2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}], (3/2) TensorTableau[{{1, 3}, {2}}]},
  TestID -> "orthogonalization" ];
VerificationTest[
  ToTensor[(TensorTableau[{{1, 2}, {3}}] + TensorTableau[{{1, 3}, {2}}])/Sqrt[6]],
  ToTensor[TensorTableau[{{1, 2}, {3}}]/Sqrt[6] + TensorTableau[{{1, 3}, {2}}]/Sqrt[6]],
  TestID -> "scalar-over-sum-distributes" ];
VerificationTest[ Head[TableauForm[(TensorTableau[{{1, 2}, {3}}] + TensorTableau[{{1, 3}, {2}}])/Sqrt[6]]] =!= TableauForm,
  True, TestID -> "tableauform-evaluates" ];
(* regression: inner product of a normalized state with itself is 1 (TensorDot must
   distribute the scalar-times-sum that ToTensor produces for a normalized tableau) *)
VerificationTest[
  With[{na = TableauNormalization[TensorTableau[{{1, 2}, {2}}]]}, TableauDot[na, na]],
  1, TestID -> "normalized-self-dot-is-1" ];

(* regression: TensorNorm must agree with the inner product even when coincident-Psi terms
   carry non-combining (symbolic) coefficients. Summing |coeff|^2 term-by-term dropped the
   cross term, giving Sqrt[2] instead of 2 for (1+x)Psi[1] at x->1. *)
VerificationTest[ TensorNorm[(1 + x) Psi[1]] /. x -> 1, 2, TestID -> "tensornorm-noncombining-coeff" ];

(* regression: inner products involving the literal 0 -- e.g. the residual of orthogonalizing
   two parallel states -- must be 0, not a leaked private iDotPsi/iToTensor symbol. *)
VerificationTest[ TensorDot[0, 0], 0, TestID -> "tensordot-zero-no-leak" ];
VerificationTest[ TableauDot[0, 0], 0, TestID -> "tableaudot-zero-no-leak" ];
VerificationTest[
  With[{d = Last @ TableauOrthogonalization[TensorTableau[{{1, 2}, {3}}], 3 TensorTableau[{{1, 2}, {3}}]]},
    {d, TableauDot[d, d]}],
  {0, 0}, TestID -> "orthogonalization-parallel-residual-zero" ];
