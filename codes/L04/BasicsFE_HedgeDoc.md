
# Basics of Finite Element Analysis

# Lecture 4
# 23.9.2026

## Warm-up and recap from the previous lecture

Welcome!

We will use this HedgeDoc during the lecture for short questions and exercises.

HedgeDoc for this lecture:  
https://pad.nixnet.services/q8-wXzI0RcKLJlKJgsLYsQ

**For the voting questions: add one `o` after the answer you think is correct.**

**Please stay in View mode whenever you are not answering a question.**  
When I ask a question, switch briefly to Edit or Both mode, add your `o`, and then return to View mode.

This helps to keep HedgeDoc responsive when many students are connected at the same time.

Please do not delete or modify other students' answers.

---

# Warm-up

## 1. How are you doing today?

Great:oooo

Good: ooooooo

A bit tired: oo

Still waking up: ooo
b
---

## 2. What did you have for breakfast?

Coffee or tea:oooo

Bread / sandwich:ooo

Porridge:oooo

Something else:o

No breakfast: oooooooo

---d

# Recap from the previous lecture

## 3. How well do you remember the previous lecture?

Very well:o

Mostly well:o

A little: oooooooooooo

Not much:o

-c--

## 4. Shape functions

For a two-node linear rod/truss element, the shape functions are

A) constant functions:

B) linear functions: oooooooooooooooo

C) quadratic functions:

D) trigonometric functions:

Not sure:

---

## 5. Shape-function property

At node 1, the shape function $N_1$

A) has the value 0:

B) has the value 1:oooooooo

C) has the value $L$:

D) depends on Young's modulus: 

Not sure: oooooo

---

## 6. Two-node rod element

The strain field of a two-node linear rod/truss element is

A) constant:oo oooo

B) linear:oo

C) quadratic:

D) cubic:

Not sure: oooooo

---

## 7. Coordinate transformation

Why do we transform the stiffness matrix of an inclined rod element from local to global coordinates?

A) Because Young's modulus changes with orientation:

B) Because elements may have different orientations in the global coordinate system: oooooooooooooo

C) Because the local stiffness matrix is incorrect:

D) Because forces can only be defined in local coordinates:

Not sure:o

---

## 8. Transformation of the stiffness matrix

Which expression is used to transform the local elemental stiffness matrix into global coordinates?

A) $\boldsymbol{K}_G = \boldsymbol{T}\overline{\boldsymbol{K}}$:

B) $\boldsymbol{K}_G = \boldsymbol{T}^{T}\overline{\boldsymbol{K}}\boldsymbol{T}$:ooooooooooo

C) $\boldsymbol{K}_G = \overline{\boldsymbol{K}}\boldsymbol{T}$:

D) $\boldsymbol{K}_G = \boldsymbol{T}^{-1}\overline{\boldsymbol{K}}$:

Not sure:oo

---

## 9. Assembly

Two rod elements share the same node.

What happens to their stiffness contributions associated with the common global degrees of freedom?

A) Only the stiffer element is used:

B) The stiffness contributions are added together: oooooooooooo

C) The stiffness contributions are multiplied:

D) One of the contributions is removed:

Not sure:o

---

## 10. Degrees of freedom

The three-rod truss example in Bonus Task 2 has four nodes and two global displacement DOFs at each node.

How many global nodal displacement DOFs does the complete system have **before applying boundary conditions**?

A) 2:

B) 4:

C) 6:

D) 8: oooooooooooo

Not sure:

---

## 11. Boundary conditions

In the three-rod truss example in Bonus Task 2, only $u_2$ and $v_2$ remain unknown after applying the prescribed zero-displacement boundary conditions.

What is the size of the reduced stiffness matrix?

A) $1 \times 1$:

B) $2 \times 2$: ooooooo

C) $4 \times 4$: o

D) $8 \times 8$:

Not sure: ooooo

---

## 12. After solving the global nodal displacements

How do we obtain the local elemental nodal displacements of element $i$?

A) $\overline{\boldsymbol{u}}^{(i)} = \boldsymbol{T}^{(i)}\boldsymbol{u}_G^{(i)}$:ooooooo

B) $\overline{\boldsymbol{u}}^{(i)} = \overline{\boldsymbol{K}}^{(i)}\boldsymbol{u}_G^{(i)}$:

C) $\overline{\boldsymbol{u}}^{(i)} = \boldsymbol{T}^{(i)T}\boldsymbol{f}_G^{(i)}$:

D) They cannot be calculated after assembly:

Not sure: ooo

---

## 13. Member forces

After the local elemental nodal displacements are known, the local elemental member forces are obtained from

A) $\overline{\boldsymbol{f}}^{(i)} = \overline{\boldsymbol{K}}^{(i)} \overline{\boldsymbol{u}}^{(i)}$:ooooooooo

B) $\overline{\boldsymbol{f}}^{(i)} = \boldsymbol{T}^{(i)} \overline{\boldsymbol{u}}^{(i)}$:

C) $\overline{\boldsymbol{f}}^{(i)} = E\overline{\boldsymbol{u}}^{(i)}$:

D) $\overline{\boldsymbol{f}}^{(i)} = \boldsymbol{K}_G^{-1} \overline{\boldsymbol{u}}^{(i)}$:

Not sure:ooo

---

# One more question before we continue

## What was the most difficult topic in the previous lecture?

Element derivation:

Shape functions:o

Coordinate transformations:o

Assembly of the global stiffness matrix: ooooooooo

Boundary conditions:o

Solving nodal displacements:

Computing member forces:

Everything was reasonably clear:

---

# Let's continue here around 8:35...


# Questions during the lecture

You can write questions here at any time during the lecture.

## Q1 I do not understand this and that...could you repeat?

## Q2 How do we know how many boundary conditions are required for different systems

## Q3

## Q4 

---

# Beam elements – quick questions

**Add one `o` after the answer you think is correct.**

---

## 14. When is a beam model a reasonable approximation?

A) When all three dimensions of the structure are approximately equal:

B) When the length is much larger than the cross-sectional dimensions:ooooooooooooo

C) Only when the beam is made of steel:

D) Only when the beam is straight:

Not sure:

---

## 15. Euler–Bernoulli beam theory

Which statement is an assumption of Euler–Bernoulli beam theory?

A) Cross-sections remain plane and perpendicular to the neutral axis:oooooooooo0o0o

B) Shear deformation is the dominant deformation mode:

C) Cross-sections deform freely during bending:

D) The beam cannot have rotations:

Not sure:

---

## 16. Euler–Bernoulli beam equation

For constant bending stiffness $EI$, which equation describes the transverse displacement $v(x)$ under a distributed load $q(x)$?

A) $EI\frac{dv}{dx}=q(x)$:

B) $EI\frac{d^2v}{dx^2}=q(x)$:

C) $EI\frac{d^3v}{dx^3}=q(x)$:oooo

D) $EI\frac{d^4v}{dx^4}=q(x)$:oooooo

Not sure:ooooo

---

## 17. Point load at the free end

Consider a cantilever beam with constant bending stiffness $EI$ and a point load $P$ at the free end.

There is no distributed load along the beam, so $q(x)=0$.

Which statement is correct?

A) The beam equation becomes $EI\frac{d^4v}{dx^4}=P$:

B) The beam equation becomes $EI\frac{d^4v}{dx^4}=0$, and the point load $P$ enters through a boundary condition at $x=L$:ooooooo

C) The beam equation becomes $EI\frac{d^2v}{dx^2}=P$: 

D) Since $q(x)=0$, the beam has no deflection:

Not sure:oooo

---

## 18. Two-node Euler–Bernoulli beam element

How many degrees of freedom does the 2D four-DOF Euler–Bernoulli beam element have at each node?

A) One – transverse displacement only:

B) Two – transverse displacement and rotation:oooooooooo
oo
C) Three – axial displacement, transverse displacement and rotation:

D) Four:

Not sure:o

---

## 19. Interpolation of the beam deflection

Why is a cubic polynomial used for the transverse displacement $v_h(x)$ of the two-node four-DOF (and six-DOF) Euler–Bernoulli beam element?

A) Because there are four nodal degrees of freedom involved in the transverse displacement interpolation: $\overline{v}_1$, $\overline{\theta}_1$, $\overline{v}_2$, and $\overline{\theta}_2$:o oo

B) Because all distributed loads are cubic:

C) Because the bending moment is always cubic:

D) Because Young's modulus varies cubically:

Not sure:ooooo

---

## 20. Rotation and deflection

How is the rotation $\theta$ approximated in the Euler–Bernoulli beam element?

A) $\theta \approx v_h$:

B) $\theta \approx \frac{dv_h}{dx}$:ooo

C) $\theta \approx \frac{d^2v_h}{dx^2}$:

D) $\theta \approx \frac{d^3v_h}{dx^3}$:

Not sure:ooooo

---

## 21. Curvature

Which quantity is related to the curvature of the Euler–Bernoulli beam?

A) $v_h$:

B) $\frac{dv_h}{dx}$:

C) $\frac{d^2v_h}{dx^2}$:o

D) $\frac{d^4v_h}{dx^4}$:

Not sure:ooooooo

---


## 23. Coupling in the 6-DOF Euler–Bernoulli beam element

Consider the following local stiffness matrix of the 6-DOF Euler–Bernoulli beam element:

$$
\overline{\boldsymbol{K}} =
\begin{bmatrix}
\frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\
0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\
-\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\
0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L}
\end{bmatrix}
$$

with the vector of elemental nodal coordinates

$$
\overline{\boldsymbol{u}} =
\begin{bmatrix}
\overline{u}_1 &
\overline{v}_1 &
\overline{\theta}_1 &
\overline{u}_2 &
\overline{v}_2 &
\overline{\theta}_2
\end{bmatrix}^{T}.
$$

What can you conclude about the coupling between the degrees of freedom from the local stiffness matrix?

A) The axial DOFs $\overline{u}_1$ and $\overline{u}_2$ are coupled with each other, and the bending DOFs $\overline{v}$ and $\overline{\theta}$ are coupled with each other, but axial and bending DOFs are not coupled:o

B) All six DOFs are coupled with each other:

C) The axial and transverse displacement DOFs are coupled, but the rotations are independent:o

D) All six DOFs are completely independent:

Not sure:ooooo

## 24. Coupling between nodes in the global stiffness matrix – a simple truss example

Consider a structure consisting of two identical two node linear rod elements:

$$
1 \;-\; 2 \;-\; 3
$$

Element 1 connects nodes 1 and 2, and element 2 connects nodes 2 and 3.

Each node has one axial displacement DOF:

$$
\boldsymbol{u}_G =
\begin{bmatrix}
u_1 & u_2 & u_3
\end{bmatrix}^{T}.
$$

The assembled global stiffness matrix is

$$
\boldsymbol{K}_G =
\frac{EA}{L}
\begin{bmatrix}
1 & -1 & 0 \\
-1 & 2 & -1 \\
0 & -1 & 1
\end{bmatrix}.
$$

What can you conclude from the non-zero and zero off-diagonal terms of the global stiffness matrix?

A) Nodes 1 and 2 are directly coupled, nodes 2 and 3 are directly coupled, but nodes 1 and 3 are not directly coupled:oooo

B) All three nodes are directly coupled with each other:

C) Only nodes 1 and 3 are directly coupled:

D) The off-diagonal terms do not contain information about the connectivity of the structure:

Not sure:oooo





## Bonus question

For a cantilever beam with a point load $P$ at the free end, the exact Euler–Bernoulli displacement solution is a cubic polynomial in $x$.

A standard two-node Euler–Bernoulli beam element also uses cubic interpolation for $v_h(x)$.

What does this suggest?

A) One beam element can represent this particular displacement field exactly:o

B) At least four beam elements are always required:

C) The finite element solution must be linear: 

D) The beam element cannot represent the analytical solution:

Not sure:oooo

---

## Should we use HedgeDoc also in a class?

A) Yes: ooo
B) No
C) Does not matter. o0oooo0

## Are contact exercises useful?

A) Yes:oo
B) No
C) Does not matter.

# End-of-lecture feedback

# The course feedback including exercises

## What is still unclear?

Write a short comment below:

- Where are bonus task 2 MatLab files? I could not ifnd them in the bonus task 2 folder in github
<span style="color:green">/bonustasks/BonusTask2/B02_Simple_truss_3rods
</span> 
- Where are the description on what to do in Exercise 1 and submission box for Bonus Task 2
<span style="color:green">The description and instructions for Exercise 1 will be provided after the next exercise session. The submission page for Bonus Task 2 is already available on Moodle. </span> 
  


