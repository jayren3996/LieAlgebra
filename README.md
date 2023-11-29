# ClassicalLieAlgebra

A **simple Lie algebra** is a Lie algebra that is non-abelian and contains no nonzero proper ideals. A finite-dimensional simple complex Lie algebra is isomorphic to either one of the following **classical Lie algebras**: 

- $A_{n-1}$: generators of $\mathrm{SU}(n)$ group,
- $B_n$: generators of $\mathrm{SO}(2n+1)$ group,
- $C_n$: generators of $\mathrm{SO}(2n)$ group,
- $D_n$: Generators of $\mathrm{USp}(2n)$ group,

or one of the five exceptional Lie algebras: $\mathrm{G}_2$, $\mathrm{F}_4$, $\mathrm{E}_6$, $\mathrm{E}_7$, and $\mathrm{E}_8$. This `Mathematica` package helps construct the linear irreducible representations of the classical Lie algebras. 

## Usage

Place the file `ClassicalLieAlgebra.wl` in a folder, create a new notebook file with extension `.nb` in the same folder, and input the import command on the first line:

```mathematica
Import[NotebookDirectory[]<>"ClassicalLieAlgebra.wl"];
```

For more examples of usage, see `Demo.nb` notebook.

## Example of su(3)

Here we discuss the SU(3) group, which has no essential difference compared to higher SU(N) groups. 

The first function of this package is to give us the specific information of the algebra generators. For example, the command to obtain the standard generators of SU(3) is:

```mathematica
MatrixForm/@Generators[SU[3]]
```

The result is 8 Gell-Mann matrices:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O1.png)

We can also obtain the **Cartan-Weyl basis** of this algebra using the following command:

```mathematica
{h, e, f} = CartanWeyl[SU[3]];
Print["H = ", MatrixForm /@ h, ", E = ", MatrixForm/@e, ", F = ", MatrixForm /@ f];
```

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O2.png)

For building representations, the most useful generator is the **Chevalley basis**, whose specific form is:

```mathematica
{h, e, f} = Chevalley[SU[3]];
Print["H = ", MatrixForm /@ h, ", E = ", MatrixForm /@ e, ", F = ", MatrixForm /@ f];
```

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O3.png)

We see that under the Chevalley basis, the generators form two coupled SU(2) groups, whose relationship can be visually expressed as a 3-level system transition:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/P1.png)

We see that the three types of Chevalley basis generators have the following effects on this 3-level system:

- $E_1, E_2$ are operators that raise the energy levels, while $F_1, F_2$ are operators that lower the energy levels.
- $E_1, F_1, H_1$ form the SU(2) algebra between $\left| 1 \right \rangle$ and $\left| 2 \right \rangle$, while $E_2, F_2, H_2$ form the SU(2) algebra between $\left| 2 \right \rangle$ and $\left| 3 \right \rangle$.
- $H_1, H_2$ correspond to the two "magnetic quantum numbers" of this 3-level system, and each energy level has an eigenvalue determined by $H_1, H_2$.
- From the eigenvalues of $H_1, H_2$, we can see that $E_1$ increases the magnetic quantum numbers by $(2,-1)$, while $E_2$ increases the magnetic quantum numbers by $(-1,2)$. Correspondingly, the operators $F_1, F_2$ decrease the corresponding magnetic quantum numbers.

This 3-level system provides the fundamental representation of SU(3). To obtain larger representations, we only need to consider the representations given by $N$ such 3-level systems under the action of generators. For an N-body system, the total generators are:
$$
H = \sum_i H_i ,\quad  E = \sum_i E_i ,\quad  F = \sum_i F_i.
$$
The total magnetic quantum numbers can be used to label different states in the representation. This is also known as the weight of the state, and states with different weight values must be orthogonal. But sometimes certain weight values in space are expressed as several linearly independent states. These states may not be orthogonal, and this is when we encounter the issue of degenerate weights. At this point, orthogonalization is required. The procedure for orthogonalization is similar to that in quantum mechanics.

### Young Tableau and Wave Function

To obtain an irreducible representation space, we actually only need to start from any wave function in this space and continuously act on the group generators to generate a closed space, which is the irreducible representation space. Therefore, the most important thing for the SU(3) irreducible representation is to find a state for each irreducible representation space.

The Young tableaux provide such a state, which is the highest-weight state of the irreducible representation. The SU(3) group irreducible representation can be represented by a Young tableau $[\lambda_1, \lambda_2]$ with no more than two rows, and the corresponding representation is denoted as $(\lambda_1-\lambda_2,\lambda_3)$. The index of the representation represents the weight value of the highest weight state or the "magnetic quantum number" of the highest weight state (determined by the eigenvalue of the Chevalley basis generator $H_i$). The tensor Young tableau of the highest weight state of each representation is to fill all the 1 in the first row of the Young diagram and all the 2 in the second row. For example, for the $(1,1)$ representation, its tensor Young tableau is:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y1.png)

The corresponding wave function of this tensor Young tableau can be implemented with the following command:

```mathematica
ct = Tableau[{{1, 2}, {3}}];
v = Psi[1, 1, 2];
TableauPermute[ct, v]
```

The output result is

```mathematica
2 Psi[1,1,2] - Psi[1,2,1] - Psi[2,1,1]
```

That is to say, for the $(l_1,l_2)$ representation of SU(3), the tensor Young tableau with the shape $[l_1+l_2,l_2]$ filled according to the above rules is the highest weight state of this representation.

Next, as long as we continuously act on the highest weight state with the generators until we cannot obtain more linearly independent wave functions, we can obtain the representation space. From here on, we can actually complete the next calculation entirely on the wave function. However, we still keep the form of the tensor Young tableau because it provides a "compact" form for the many-body wave function. At the same time, the tensor Young tableau intuitively displays some permutation symmetries of the wave function. For example, since the elements of the Young operator are all anti-symmetric, the same number cannot appear in the same column of the tensor Young tableau, otherwise the wave function is 0.

### The (1, 1) representation of su(3)

Now, we start from the highest weight state and act on two generators under the Chevalley basis on this wave function. The new wave function obtained by the action represents a "transition" process under the action of the generators, and we can use lines to represent a "transition" process (a single line represents the action of $F_1$, and a double line represents the action of $F_2$). We temporarily do not discuss the size of the transition matrix elements here. The results of the first-level "transition" are:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/T1.png)

We labeled the weight of each state next to each Young tableau (the weight calculation can be obtained by adding up the magnetic quantum numbers of each small cell or by using the lowering operator to obtain the weight change of the weight).

Now let's consider level 2 and pay attention to the state on the left:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y4.png)

Under the action of $F_1$, it is annihilated (the same number appears in the same column), and under the action of $F_2$, it becomes a superposition state of two tensor Young diagrams:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y5.png)

On the right-hand side, the state is annihilated under the action of $F_2$ (no digit 2 appears in the tensor Young diagram), and the state obtained under the action of $F_1$ is:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y7.png)

The second table is not a regular Young tableau, but can be symmetric using the symmetry of the Young operator (1 exchanges with 2,3).

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y8.png)

Converted to regular tensor Young diagram, finally:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y9.png)

At this point, the weights of the two wave functions obtained by level 3 are $(0,0)$, which indicates a situation of heavy weight. At this time, we need to orthogonalize the weight space. Just like the orthogonalization in quantum mechanics, the choice of the basis is not unique. The customary selection is to retain the state generated by the generator with smaller index, and orthogonalize the other states relative to this state. Here, it corresponds to retaining the state on the right. Now, we can orthogonalize these two Young diagrams into wave functions. This mechanical procedure can be completed using the Mathematica package. First, establish two tensor Young diagrams:

```mathematica
a = 2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}];
b = TensorTableau[{{1, 3}, {2}}] + TensorTableau[{{1, 2}, {3}}];
```

The orthogonalization function is:

```mathematica
{c, d} = TableauOrthogonalization[a, b];
```

Finally, we can print the results:

```mathematica
TableauForm /@ {c, d}
```

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O4.png)

In this way, we have drawn the structure diagram of the (1,1) representation of the SU(3) Lie algebra, and we see that this is an 8-dimensional space, and 8 nodes give a set of orthonormal bases for this space. At the same time, the transition behavior between the generating elements in this space basis state may also be determined by this structural diagram.

Now, there is one more question, which is the size of the transition matrix element between the generating elements among these nodes. This is actually similar to a quantum mechanics problem. The coefficient of the transition is largely determined by the normalization of each state. Here we analyze a specific transition process, which is the one in the above figure:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/T4.png)

We start from the state:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/Y4.png)

and denote it as $\left| a \right\rangle$, first normalize it by the following command:

```mathematica
a = TensorTableau[{{1, 2}, {2}}];
na = TableauNormalization[a];
Print["|a> = ", TableauForm[na]];
```

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O5.png)

The result of the action of the generating element $F_2$ is:

![img](https://raw.github.com/jayren3996/LieAlgebra/master/pics/O6.png)

We denote this state as $\left| b \right\rangle$, and the two basis states $\left| c \right \rangle,\left| d \right \rangle$​. First, we normalize the states $\left| c \right \rangle$ and $\left| d \right \rangle$ and calculate their inner product separately. This calculation can be performed using the following command (which also displays the normalized states):

```mathematica
b = TensorTableau[{{1, 2}, {3}}]/Sqrt[6] + TensorTableau[{{1, 3}, {2}}]/Sqrt[6];
c = TensorTableau[{{1, 3}, {2}}];
d = 2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}];
nc = TableauNormalization[c];
nd = TableauNormalization[d];
Print["|c> = ", TableauForm[nc], ", |d> = ", TableauForm[nd], ", |d> = ", TableauForm[nd],", <d|b> = ", TableauDot[nd, b]];
```

The corresponding wave function inner product is the transition matrix element.

