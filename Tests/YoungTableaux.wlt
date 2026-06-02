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
