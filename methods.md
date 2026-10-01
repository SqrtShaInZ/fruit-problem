## Methods used

### Preliminary

A PARI/GP script was used for $n\le20000$. Later an equivalent Sagemath script was used for $20001\le n\le 100000$. They all feature `ellrank` with `effort=3`.

Approximately $65\%$ entries of the entire database were able to be filled completely after this round. The missing generators of the rest entries are denoted as `[]`. Lines representing those entries could look like:

```
615 [[], []]
616 [[]]
617 []
618 [[-276/5041, 21215844/357911], []]
```

Note that there might be extra "missing generators" due to $\mathrm{Sha}[2]$, because `ellrank` only computes the $2$-Selmer rank. To fix this problem, Magma's `ThreeDescent` rank bounding was used in the Sagemath script, together with `analytic_rank_upper_bound` with $\Delta=2.5$.

### $2$-Selmer Rank $\ge 4$

Since very few curves in the whole database has rank $\ge 3$, all of them are processed manually (though lots of them already have no missing generators).

There was a critical bug in the Sagemath script that caused many fake rank $\ge3$ curves to emerge due to inaccurate translation from PARI/GP. The bug was later fixed and all such indices have been rerun.

### $2$-Selmer Rank $3$ with one missing generator

According to Cassels-Tate pairing, $\mathrm{Sha}$'s contribution to the $2$-Selmer rank must be an even number, hence the algebraic rank is $2$. One can verify in another way through root number, assuming the parity conjecture.

`ComputeGeneratorFull` routine was used. Before that, the `best` function in PARI/GP script was used to determine the best curve based on the known generator and BSD formula. Note that `best` is only designed for Rank $2$

#### How `best` works

`best` calculates the BSD factor `c` in the BSD formula

$$
\frac{L^{(2)}(E,1)}2=c(E)\times\mathrm{Reg}(E)\times\mathrm{Sha}(E)
$$

and determines which curve(s) in the isogeny class has smallest $\mathrm{Reg}(E)\times\mathrm{Sha}(E)$ as predicted by the formula, using the fact that the L-function is the same across an isogeny class.

* When there is only one curve with smallest $\mathrm{Reg}(E)\times\mathrm{Sha}(E)$, it performs `ellrank` with `effort=3` on it so the known generator can be found.

* When there are two such curves, `ellrank` with `effort=3` is used to find generators of both curves, and `best` automatically chooses the one with larger height to guarantee that the second generator has smaller height (contributions from Tate-Shafarevich groups are completely ignored here). Sometimes only one of the two `ellrank` calls returns a point, so we have to use `ellisogeny` formulas to find the first generator for the other curve.

* When all curves in the isogeny class has the same $c$, the routine is rather stupid: It emulates the database generation by calling `findgen_worker` with `effort=-600000` (which corresponds to `ellrank`'s `effort=3` and avoids analytic calculations), uses `ellisogeny` to find a point on each curve in the isogeny class, saturates them, and returns the point with largest height. It is true that more efficient code can be introduced, but the current version is working well.

#### How `ComputeGeneratorFull` works

It is an adaptive algorithm like `MordellWeilShaInformation`. Specifically, it follows the order how the latter executes the descent formulas, but with two major changes:

* 6- and 12-descent are implemented, and they are used after 3-descent calls. This makes finding higher generators more accessible, particularly those out of reach of 8-descent.

* 8-descent is not performed on 3-isogenous or 6-isogenous curves. This is because in those cases 6-descent and 12-descent are way cheaper than 8-descent, which involves highly sophisticated algorithms to find solution to conics.

It is recursive, which means it appends the newly found generator to `known_gens` parameter, then calls itself. The recursion boundary is `#known_gens==rank`, when all the generators are known, which then triggers saturation.

### $2$-Selmer rank $3$ with two missing generators

Suppose the entry is `n`.

In those cases the algebraic rank might be $0$. So a pass of `TwoPowerIsogenyDescentRankBound(mkc(n))` is computed, bounding some curves to rank $0$. Another pass consists of Sagemath's `mkc(n).analytic_rank_upper_bound(adaptive=True, max_Delta=3.0)`

When the above function converges to a value around $2.2$, then the isogeny class is expected to be rank $2$, and we jump to the "Knowing the curve has rank $2$" section.

For those entries that are still like `[[], []]`, this following snippet is executed:

```
Crv3 := ThreeDescent(mkc3(n));
[#pIsogenyDescent(Crv3[i], mkc(n), mkc3(n)) : i in [1..#Crv3]];
```

When the snippet successfully executes, and the vector is $[0,0,0,0]$, the rank is definitely $0$, since none of the $3$-covers lift to a $9$-cover. This is cheaper than 8-descent.

There is one final pass that might prove the rank is $0$; it relies on the following lemma.

> **Lemma.** $E(\mathbb{Q})$ is an elliptic curve with $E(\mathbb{Q})[2]\cong \mathbb{Z/2Z}$ and $2$-Selmer rank $3$. Let $\pi_2 : \mathrm{Sel}^{(2)}(E)\to\mathrm{Sel}^{(2)}(E)/E(\mathbb{Q})[2]$ be the quotient map, and the nonzero elements of this quotient are $\{HE_1,HE_2,HE_3\}$. The $4$-coverings $\{C4_{i,1},C4_{i,2}\}\sub \mathrm{Sel}^{(4)}(E)$ map to $HE_{i}$ under the projection $\pi : \mathrm{Sel}^{(4)}(E)\to\mathrm{Sel}^{(2)}(E)/E(\mathbb{Q})[2]$. If $C4_{i,1},C4_{i,2}$ are all elements of $\mathrm{Sha}(E)[4]$ for some $i\in[1,3]$, then $\mathrm{rank}(E(\mathbb{Q}))=0$.
> 
> **Proof.** From the exact sequence $0 \to E(\mathbb{Q})/2E(\mathbb{Q}) \to \text{Sel}^{(2)}(E) \to \text{Sha}(E)[2] \to 0$, we have:
> 
> $$
> \dim_{\mathbb{F}_2} \text{Sel}^{(2)}(E) = \text{rank}(E(\mathbb{Q})) + \dim_{\mathbb{F}_2} E(\mathbb{Q})[2] + \dim_{\mathbb{F}_2} \text{Sha}(E)[2]
> $$
> 
> Substituting $\dim_{\mathbb{F}_2} \text{Sel}^{(2)}(E) = 3$ and $\dim_{\mathbb{F}_2} E(\mathbb{Q})[2] = 1$:
> 
> $$
> 3 = \text{rank}(E(\mathbb{Q})) + 1 + \dim_{\mathbb{F}_2} \text{Sha}(E)[2] \\
\implies \text{rank}(E(\mathbb{Q})) + \dim_{\mathbb{F}_2} \text{Sha}(E)[2] = 2
> $$
> 
> Since $C4_{i,1}$ and $C4_{i,2}$ represent non-trivial elements in $\text{Sha}(E)[4]$, their projection $HE_i$ generate a subspace of $\text{Sha}(E)[2]$. Therefore $\dim_{\mathbb{F}_2} \text{Sha}(E)[2] \ge 1$. The parity implication of the Cassels–Tate theorem immediately yields $\dim_{\mathbb{F}_2} \text{Sha}(E)[2]\ge 2$.
> 
> It follows that $\text{rank}(E(\mathbb{Q})) \le 2 - 2 = 0$. $\square$

The procedure itself is straightforward. Proving that $4$-covers are elements of $\mathrm{Sha}(E)[4]$ can be achieved through `EightDescent`.

#### Knowing the curve has rank $2$

Compute `exp1(n)`, `exp2(n)`, `exp3(n)`, `exp6(n)` (those are PARI/GP functions). Choose the one(s) with smallest `exp` and proceed. Since `exp` is designed as $\frac 1{\texttt{ellbsd}(E)}$, a smaller `exp` means a smaller regulator, assuming $\mathrm{Sha}$ holds invariant.

Call `MordellWeilShaInformation` on each of the chosen curves.

* When one call returns two generators, the process terminates; the full Mordell-Weil basis has been found.

* When one call returns a generator, stop further `MordellWeilShaInformation` calls and use algorithms inside the "$2$-Selmer Rank $2$ with one missing generator" section.

* When none of the calls return a generator, run `ellanalyticrank`.
  
  * If the rank is confirmed to be $0$, process halts.
  
  * If the rank is confirmed to be $2$, continue with higher descent methods. `ellanalyticrank` also provides a clear view of the regulator, combined with known information of $\mathrm{Sha}$ from $4$-descent.

### $2$-Selmer rank $2$ with no known generators

First, the best curve is determined using BSD formula again (there is a dedicated function `BSDTermsEasyQ`). $4$-descents are performed before a `ConjecturalRegulator` call. This preliminary $4$-descent acts as a filter that can solve around $25\%$ of the cases.

According to the regulator, several descent methods are performed, in the following order:

* When regulator $\le135$, $4$-descents are performed.

* When regulator $\le540$, and `NoEightDesc` is false, $8$-descents are performed. It can fall back to $12$-descent.

* When regulator $\le 240$, $6$-descents are performed.

* $12$-descents are performed.

For the automated phase, the regulator is limited to $1500$. This prevents high-regulator curves stalling the automated pipeline.

#### Estimating $\mathrm{Sha}$ from descent methods

The BSD-predicted regulator requires $|\mathrm{Sha}(E)|$. Its $2$-part is bounded below using the number $N$ of $4$-coverings returned by `FourDescent(E)`.
Assume $E(\mathbb{Q})[2]\cong\mathbb{Z}/2\mathbb{Z}$ and $\mathrm{rank}\,E(\mathbb{Q})=1$. Write

$$
N = 2^p\,(2^q-1)
$$

Then $p=\dim\mathrm{Sha}(E)[2]$, $q$ is the dimension of the image of $\mathrm{Sel}^{(4)}(E)$ in $\mathrm{Sel}^{(2)}(E)/E(\mathbb{Q})[2]$ (the kernel of the Cassels–Tate pairing), and

$$
|\mathrm{Sha}(E)[4]| = 2^{p+q-1}
$$

The factorisation is unique, and two consistency checks come for free: $q-1\le p$
(since $2\,\mathrm{Sha}[4]\subseteq\mathrm{Sha}[2]$) and $q-1\equiv p \pmod 2$ (the pairing is alternating). For example $N=28$ is $(p,q)=(2,3)$, which corresponds to $\mathrm{Sha}[4]\cong(\mathbb{Z}/4\mathbb{Z})^2$.

The $3$-part is estimated from the number of $3$-coverings (`MyThreeDescent`),
as before.

#### Faster $3$-descent for Rank $1$ elliptic curves

To accelerate 3-descent on curves admitting 3-isogenies, a custom pipeline (`MyThreeDescent`) replaces the standard Magma `ThreeDescent` routine. Standard 3-descent can stall when computing full 3-coverings directly; exploiting 3-isogenous pairs $E \xrightarrow{\phi} E_\phi$ significantly lowers the computational algebra overhead.

- For curves $E$ with a 3-isogeny to $E_\phi$, the algorithm first computes `ThreeIsogenyDescent(E)`, yielding candidate coverings and maps.

- Maps are composed with the dual isogeny $\hat{\phi}$.

- For each 3-covering $C_{\phi, i}$ of the target curve $E_\phi$, it pulls the covering back to $E$ using `pIsogenyDescent(Crv3_phi[i], E, E_phi)`.

### Saturation

This is accomplished using PARI/GP's `ellsaturation` with bound $1000$. The saturated basis is then checked against BSD-predicted regulator (for Rank $1$ curves).

### Rooms for improvement

* Several steps are conditional, and rely on one or more of the following conjectures: the parity conjecture, Generalized Riemann Hypothesis, and the Birch and Swinnerton-Dyer formula. Specifically, $\mathrm{Sha}$ is assumed to be finite for every elliptic curve.

* The heuristic can definitely misidentify the optimal curve for descents, especially those with rare $\mathrm{Sha}$ structures. But taking care of those $\mathrm{Sha}$ is too complicated, and is omitted in my script.

* The `pIsogenyDescent` snippet can still throw errors because sometimes `ThreeDescent` returns full 3-covers that don't possess a rational flex.
