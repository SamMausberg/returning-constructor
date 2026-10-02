import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
src = root / "source"
aux = (root / "_build/paper.aux").read_text()
labels = dict(re.findall(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}", aux))
items = [
    (
        "thm_entropy",
        "Closed finite memory, independent input state $\\tau^{\\otimes n}$, $D\\ge2$, $m=\\log_2D$, and endpoint return error $0\\le\\delta\\le1/2$. Initial memory may be mixed. The initial state is a tensor product of inputs and memory.",
        r"Audenaert~\cite{audenaert}, Theorem~1 (entropy continuity). The closed-unitary entropy balance is part of the statement to be checked.",
    ),
    (
        "cor_exactreturn",
        "The same closed model with independent $\\tau^{\\otimes n}$ inputs, arbitrary finite $D\\ge1$, and exact endpoint return $\\sigma_n=\\sigma_0$.",
        r"The preceding entropy statement. No additional source-dependent assumption.",
    ),
    (
        "thm_forwardoptimal",
        "$n\\ge1$ pure-zero inputs, $m\\ge0$ memory qubits, arbitrary mixed initializer independent of the inputs, and one fixed sequential unitary. The target joint state is $\\tau^{\\otimes n}$.",
        r"Boes et al.~\cite{boes}, Lemma~1 and Sec.~III.A (dephasing mechanism); Lie--Jeong~\cite{lie} (catalytic correlations); Douglas~\cite{douglas}, Proposition~7.2 (counted streaming). These sources identify antecedents of the construction.",
    ),
    (
        "thm_dyadic",
        "$n\\ge1$ independent maximally mixed inputs, an $m$-qubit memory with $m\\ge2n$, freely chosen mixed initializer, and fixed cyclic-swap unitary. The cost uses return after every prefix.",
        r"Ng et al.~\cite{ng}, Eqs.~(2)--(3) (endpoint catalyst). The prefix-return equality is a statement of the present paper.",
    ),
    (
        "prop_spectralreturn",
        "Independent $\\tau^{\\otimes n}$ inputs, arbitrary initial $m$-qubit memory, exact joint output $p^{\\otimes n}$, and endpoint return tolerance $\\delta>0$. A global closed unitary is allowed, so the bound also covers sequential devices.",
        r"Ng et al.~\cite{ng}, Eqs.~(2)--(3) (endpoint catalyst optimum).",
    ),
    (
        "thm_moments",
        "Forward coherent homogenizer: $p^{\\otimes n}$ inputs, initial reservoir $\\tau^{\\otimes M}$, $M\\ge1$, $0<\\eta\\le\\pi/2$, and $u=\\sin^2\\eta$. This statement uses only the first visit.",
        r"Ziman et al.~\cite{ziman}, partial-swap homogenizer definition; PRL~\cite{prl}, homogenizer model. The moment identities and constants are statements of the present paper.",
    ),
    (
        "cor_singlecost",
        "The same forward homogenizer at one visit; coupling may be optimized. The two tolerances satisfy $0<\\epsilon,\\delta\\le1/4$, and reservoir size is an integer.",
        r"The moment statement and the original partial-swap model~\cite{ziman}. No source supplies the displayed size constants.",
    ),
    (
        "thm_transport",
        "Arbitrary excitation-conserving closed unitary, initial state $p^{\\otimes n}\\otimes\\tau^{\\otimes M}$, unconditional output marginals within $\\epsilon<1/2$ of $\\tau$, and endpoint return within $\\delta>0$.",
        r"The excitation-conserving partial swap is defined in Refs.~\cite{ziman,scarani}. The bound is stated for the larger unitary class specified here.",
    ),
    (
        "lem_collective",
        "Forward homogenizer with $n,M\\ge1$, product initializer $\\tau^{\\otimes M}$, pure-zero inputs, and prescribed $0<\\eta\\le3/(8n^3)$.",
        r"Giovannetti--Palma~\cite{gp}, Eqs.~(4)--(7), and Ref.~\cite{cascade}, Eqs.~(5)--(8), for the column representation. The weak-coupling inequality is a statement of the present paper.",
    ),
    (
        "thm_manyupper",
        "The same forward homogenizer, $n\\ge1$, $0<\\epsilon,\\delta\\le1/4$, with coupling chosen according to one of the two designs. Return is required at every prefix; output accuracy is marginal.",
        r"Violaris--Marletto~\cite{vm}, Sec.~3.3, and Beever et al.~\cite{beever}, Sec.~4 (finite-horizon resource formulations). The numerical bounds are those stated here.",
    ),
    (
        "thm_singlet",
        "Forward coherent homogenizer, initial reservoir $\\tau^{\\otimes M}$, pure-zero inputs, $M\\ge1$, and fixed $0<\\eta\\le\\pi/2$. The obstruction uses at least two outputs and includes return after the first use.",
        r"The coherent partial-swap definition~\cite{ziman} and the column-channel representation~\cite{gp,cascade}. The singlet identity is a statement of the present paper.",
    ),
    (
        "lem_relaxation",
        "A finite block of $k\\ge1$ qubits, arbitrary possibly correlated block operator $X$, a fresh qubit in any state $b$ for each sweep, and fixed $0<\\eta\\le\\pi/2$.",
        r"Giovannetti--Palma~\cite{gp}, Eqs.~(4)--(7); Ref.~\cite{cascade}, Eqs.~(5)--(8); Ziman et al.~\cite{ziman}, homogenizer model. No reservoir factorization assumption.",
    ),
    (
        "thm_limits",
        "Product input and reservoir initializers, fixed coupling $0<\\eta\\le\\pi/2$, and qubit states satisfying $0<f=F(a,b)<1$. The numerator is current-output infidelity, with squared Uhlmann fidelity in both places.",
        r"Violaris--Marletto~\cite{vm}, Eqs.~(8)--(11), and thesis~\cite{thesis}, Eqs.~(5.1)--(5.4), for the ratio and orders of limits. The Letter uses a different numerator and set infimum.",
    ),
    (
        "prop_capable",
        "$M\\ge2$, strict partial swap $q=\\cos^2\\eta\\in(0,1)$, forward task $p\\to\\tau$, first-use tolerance $e_M=(1-\\sqrt{1-q^{2M}})/2$, and capability defined with squared fidelity.",
        r"PRL / arXiv:2009.14649v2~\cite{prl}, Eqs.~(3)--(7), (12)--(13). These are the definitions being tested, rather than a restricted designated-initializer condition.",
    ),
    (
        "prop_entropy-density",
        "States $a,b$ on the same finite-dimensional one-site Hilbert space, $S(a)\\ne S(b)$, independent product inputs $a^{\\otimes n}$, and fixed finite memory dimension $D$. The initial memory is independent of those inputs and may vary with the horizon.",
        r"Ryb\'{a}r--Ziman~\cite{rybar}, discussion of repeatability and finite-memory entropy constraints. The strong joint trace-error conclusion is the statement to be checked.",
    ),
]
# Include all shared definitions and the ancillary equation definitions referenced by statements.
header = r"""\documentclass[11pt]{article}
\usepackage[letterpaper,margin=1in]{geometry}
\input{preamble.tex}
\usepackage[numbers,sort&compress]{natbib}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\hypersetup{pdftitle={Statement-only referee packet: Memory return in quantum machines},pdfauthor={Samuel Mausberg}}
\begin{document}
\begin{center}
{\Large Statement-only referee packet}\\[5pt]
{\large Memory return in quantum machines}\\[5pt]
Samuel Mausberg \qquad October 2, 2026
\end{center}

This packet collects the mathematical statements and their hypotheses, without proofs, proof sketches, or numerical evidence. The source locators distinguish definitions under examination from prior mechanisms and imported inequalities. A cited source need not contain the theorem stated here. The theorem numbers agree with the accompanying paper.

\section*{Conventions}
A closed sequential device applies the same unitary to each fresh input qubit and its retained memory, then releases the qubit. The memory contains every accessible clock, buffer, random seed, and workspace. There are no resets or discarded internal registers. Initial memory is independent of the input product state. Put $p=\ket0\bra0$, $\tau=I/2$, $D=2^m$, and
\[
T(\rho,\sigma)=\tfrac12\|\rho-\sigma\|_1,
\qquad F(\rho,\sigma)=\left(\operatorname{Tr}\sqrt{\sqrt\rho\sigma\sqrt\rho}\right)^2.
\]
Entropy and $h_2$ use base-two logarithms; $\ln$ is natural. Write $\Omega_n$ for all outputs, $\rho_j$ for output $j$, and $\sigma_j$ for retained memory. Marginal accuracy is $\max_{j\le n}T(\rho_j,b)\le\epsilon$; joint accuracy is $T(\Omega_n,b^{\otimes n})\le\epsilon$. Prefix return is $\max_{j\le n}T(\sigma_j,\sigma_0)\le\delta$, whereas endpoint return checks only $j=n$. A superscript $\mathrm{ret}$ on a memory cost includes prefix return. Exact pure output marginals already fix their joint state.

The coherent homogenizer starts in $a^{\otimes n}\otimes b^{\otimes M}$. Every input visits cells $1,\ldots,M$ through $U_\eta=\cos\eta\,I+i\sin\eta\,\mathsf S$, where $\mathsf S$ swaps the two qubits. Put $u=\sin^2\eta$ and $q=1-u$. The channel $\Phi_{b,k}$ sweeps a fresh ancilla in $b$ through a $k$-qubit block and traces out that ancilla.

For the limit statements, define
\[
 R_{M,n}=\frac{1-F(\rho_{M,n},b)}{F(\sigma_{M,n},b^{\otimes M})}.
\]
The single-parameter limits referenced there are $\lim_M R_{M,n}=0$ at fixed $n$ and $\lim_n R_{M,n}=(1-f)f^{-M}$ at fixed $M$. For the capability statement, $\mathcal O_p$ is the first-output channel and
\[
V_e=\{\omega:F(\mathcal O_p(\omega),\tau)\ge1-e\},\qquad
S_e(n)=\inf_{\omega\in V_e}F(\omega,\Phi_{p,M}^{n}(\omega)).
\]

\newpage
\section*{Theorems and corollaries}
"""


def clean(s):
    s = re.sub(r"(?m)^[ \t]*\\label\{[^}]*\}[ \t]*\n", "", s)
    s = re.sub(r"\\label\{[^}]*\}", "", s)
    for cmd in ["eqref", "cref", "Cref", "ref"]:

        def conv(m):
            k = m.group(1)
            if k not in labels:
                sys.exit(f"make_packet.py: unresolved paper label {k!r}; build the paper first")
            v = labels[k]
            if cmd == "eqref":
                return "(paper Eq.~" + v + ")"
            typ = "statement" if k.startswith(("thm:", "lem:", "cor:", "prop:")) else "section"
            return "paper " + typ + "~" + v

        s = re.sub(r"\\" + cmd + r"\{([^}]+)\}", conv, s)
    return s


out = [header]
for name, h, loc in items:
    st = (src / f"statements/{name}.tex").read_text()
    out.append(
        "\n\\par\\medskip\n\\noindent\\begin{minipage}{\\linewidth}\n\\emph{Hypotheses.} "
        + h
        + "\n\n"
    )
    out.append(clean(st))
    out.append(
        "\n\\noindent\\emph{Source locators and role.} " + loc + "\n\\end{minipage}\n\\par\n"
    )
# source keys list useful indexed by contribution without proof content
out.append(r"""
\section*{Additional statements used in the paper}

\subsection*{A. Circuit exchange and reservoir-fidelity floor}
For the coherent grid with independent product initializers, all $M,n\ge1$, and fixed partial swap,
\[
 \Omega_{M,n}=\Phi_{b,n}^{M}(a^{\otimes n}),\qquad
 \sigma_{M,n}=\Phi_{a,M}^{n}(b^{\otimes M}),\qquad
 F(\sigma_{M,n},b^{\otimes M})\ge F(a,b)^n.
\]
The grid exchange holds for arbitrary two-wire gates when only disjoint-wire pairs are reordered. Source locators: Ref.~\cite{gp}, Eqs.~(4)--(7), and Ref.~\cite{cascade}, Eqs.~(5)--(8). The fidelity floor is the paper's stated consequence for the product boundary model.

\subsection*{B. Pure initial memory}
For $n$ independent maximally mixed inputs, pure initial memory of any finite dimension, and an arbitrary global closed unitary,
\[
 F(\sigma_n,\sigma_0)\le2^{-n}+e_{\rm all},\qquad
 e_{\rm all}=1-\bra{0^n}\Omega_n\ket{0^n}.
\]
Under joint error $\epsilon$, $e_{\rm all}\le\epsilon$; under marginal error $\epsilon$, $e_{\rm all}\le n\epsilon$. At one use, return requires $\delta\ge1/2-\epsilon$. For the reverse homogenizer this tradeoff is exact:
\[
 T(\rho_{M,1},p)=q^M/2,\qquad
 T(\sigma_{M,1},p^{\otimes M})=(1-q^M)/2.
\]
Source dependence: the closed-unitary and homogenizer definitions above; no stronger general-initializer conclusion is assumed.

\subsection*{C. Prescribed-coupling intervals}
For one forward homogenizer use and $0<\epsilon,\delta\le1/4$, put $\ell=\ln[1/(2\epsilon)]$. At a prescribed $u=\sin^2\eta$, necessary conditions are
\[
 \ell/[-\ln(1-u)]\le M\le128\delta^2/u^2,
\]
and sufficient conditions are $\ell/u\le M\le4\delta^2/u^2$, with integer feasibility understood. For $n$ marginal uses, $L=\ln(2/\epsilon)$, and $0<\eta\le3/(8n^3)$, the conditions
\[
 2L/\eta^2\le M\le4\delta^2/(n^2\sin^4\eta)
\]
are sufficient for prefix return and output accuracy. The one-use and transport inequalities remain necessary. Sources: the model definitions in Ref.~\cite{ziman}; finite-horizon formulations in Refs.~\cite{vm,beever}. The constants belong to this paper's statements.

\subsection*{D. The fixed-horizon collective limit}
For the forward homogenizer at fixed $n$, set $\eta=M^{-\alpha}$ with $1/3<\alpha<1/2$. Then
\[
 \Omega_{M,n}\longrightarrow\frac{\Pi_{\rm sym}}{n+1},\qquad
 \max_{j\le n}T(\sigma_{M,j},\tau^{\otimes M})\longrightarrow0.
\]
Here $\Pi_{\rm sym}$ projects onto the symmetric $n$-qubit subspace. The limiting output marginals are $\tau$, while its joint trace distance from $\tau^{\otimes n}$ is $1-(n+1)/2^n$. Source dependence: the collision-grid representation~\cite{gp,cascade}; this is a fixed-$n$ statement.

\subsection*{E. Approximate reverse outputs and prefix return}
For arbitrary mixed initial memory and $0\le\epsilon<1/2$, a device exists with
\[
 m\le1+\max\{2n,\lceil wn/\delta\rceil+n-2\},
\]
where $w=1-\epsilon$ for joint accuracy and $w=1-2\epsilon$ for marginal accuracy, for positive small $\delta$. In particular,
\[
 \lim_{\delta\downarrow0}\limsup_{n\to\infty}
 \frac{\delta m_-^{\mathrm{ret,joint}}(n,\epsilon,\delta)}n=1-\epsilon.
\]
The corresponding marginal coefficient lies between $1-h_2(\epsilon)$ and $1-2\epsilon$. Source locators: Audenaert~\cite{audenaert}, Theorem~1; Ng et al.~\cite{ng}, Eqs.~(2)--(3), for related endpoint catalysis. The statements here require prefix return.

\subsection*{F. Output-only optima}
With return unrestricted and arbitrary mixed initial memory, the optimum joint trace errors are
\[
 E_+^*(n,m)=\max\{0,1-2^{2m-n}\},\qquad
 E_-^*(n,m)=\max\{0,1-2^{m-n}\}.
\]
They also hold for squared-fidelity infidelity. With pure initialization, exact joint service in either direction requires $n$ memory qubits. Exact marginal forward service requires one qubit. At fixed $0<\epsilon<1/2$, reverse marginal service has cost
\[
 m_-^{\mathrm{marg}}(n,\epsilon)=n[1-h_2(\epsilon)]+o(n).
\]
Source locators: Ryb\'{a}r--Ziman~\cite{rybar} for repeatability, Boes et al.~\cite{boes}, Lemma~1, for catalytic dephasing, and Schulman--Vazirani~\cite{sv} for reversible compression.

For the homogenizer with coupling optimized and return unrestricted, $M=n$ attains exact joint service by full swap. Necessary marginal bounds are $M\ge n(1-2\epsilon)$ forward and $M\ge n[1-h_2(\epsilon)]$ reverse. Necessary joint bounds are
\[
 M\ge n-\sqrt{8n\ln[2/(1-\epsilon)]},\qquad
 M\ge\lceil n+\log_2(1-\epsilon)\rceil,
\]
respectively, for $0\le\epsilon<1$. At vanishing error both directions have memory rate one. Strict partial swaps have positive finite-$M$ first-use trace error $q^M/2$. Source: partial-swap definitions in Refs.~\cite{ziman,scarani}; finite-horizon resource questions in Refs.~\cite{vm,beever}.

\subsection*{G. Capability variants and limit boundaries}
For $q\ge1/2$, the initializer $(I-p^{\otimes M})/(2^M-1)$ is in $V_{e_M}$, has zero overlap with $p^{\otimes M}$, and has trace distance $2^{-M}$ from $\tau^{\otimes M}$. In the reverse direction, the capability-set infimum tends to $2^{-M}$. If the infimum is replaced by the designated product initializer while retaining the first-use numerator, the forward and reverse fixed-$M$ long-use ratios equal
\[
 2^M\frac{1-\sqrt{1-q^{2M}}}{2},\qquad \frac{(2q)^M}{2}.
\]
Their reservoir-size limits separate for $1/2<q<1/\sqrt2$. The current outputs still tend to the incoming state in both directions. For orthogonal pure endpoints and strict partial swaps, the one-use current-error ratio is one; at full swap it is undefined. Source locators: PRL~\cite{prl}, Eqs.~(3)--(7), and supplement~\cite{supp}, Eqs.~(9)--(10).

\subsection*{H. Invariant marginal capability}
For independent prescribed inputs, define $V_\epsilon=\{\sigma:T(\mathcal O_a(\sigma),b)\le\epsilon\}$ and $C_\epsilon=\bigcap_{j\ge0}\Phi_a^{-j}(V_\epsilon)$. Then $C_\epsilon$ is exactly the set of memory states serving every horizon at marginal tolerance $\epsilon$, and $\Phi_a(C_\epsilon)\subseteq C_\epsilon$. The inverse notation means set preimage. Source locators: Violaris--Marletto~\cite{vm}, Sec.~2.1; thesis~\cite{thesis}, Secs.~3.3 and 5.3. Foundational interpretation: Deutsch~\cite{deutsch}, Secs.~1.1 and 3.14, and the introductory definitions in Refs.~\cite{information,thermo}.

\section*{Version-specific source statements to check}
The PRL v1 preprint~\cite{prlv1}, p.~2, uses a supremum and one two-variable limit sign; p.~3 takes uses first in its homogenizer statements. Version~2 / the Letter~\cite{prl}, Eqs.~(3)--(7), (12)--(13), uses squared fidelity, a first-use numerator, an infimum, and uses first. The publisher supplement~\cite{supp}, Eqs.~(9)--(10), (17)--(25), substitutes the designated initializer, uses $S(n)\approx S(1)^n$, and invokes weak-coupling/product-reservoir approximations. Ref.~\cite{vm}, Eqs.~(8)--(13), uses a current-error numerator; its preprint leaves the joint limit unspecified, while the published Eq.~(11) takes reservoir size first. The thesis~\cite{thesis}, Eq.~(5.4), takes reservoir size first. Ref.~\cite{beever} treats finite-horizon marginals and individual cells, with Bloch distance $2T$. Ref.~\cite{tests}, Sec.~3.2, introduces no replacement limit theorem.

\makeatletter\immediate\write\@auxout{\string\citation{apsrev42Control}}\makeatother
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}
""")
text = "".join(out)
# Keep each ancillary statement with its conditions and source locators.
text = re.sub(
    r"(\\subsection\*\{[A-H]\..*?)(?=\\subsection\*|\\section\*|\\nocite|\\makeatletter)",
    lambda m: (
        "\\par\\medskip\n\\noindent\\begin{minipage}{\\linewidth}\n"
        + m.group(1)
        + "\n\\end{minipage}\n"
    ),
    text,
    flags=re.S,
)
(src / "referee_packet.tex").write_text(text)
