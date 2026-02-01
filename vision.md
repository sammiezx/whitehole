Below is an ambitious, high-level plan that could realistically become:

a thesis-level tool,

a publishable computational physics / GR visualization paper, or

an open research platform for causal structure and horizon formation.

I’ll give you:

The scientific vision

Core conceptual pillars

System architecture

Mathematical engine

Visualization philosophy

Research-level extensions

What makes this genuinely new

A staged execution roadmap

This is the kind of plan you would pitch to a supervisor or funding committee.

1. Scientific Vision (what problem this actually solves)

There is currently no tool that lets a researcher compute and visualize global causal structure from dynamical spacetime data.

Textbooks:

draw Penrose diagrams after the fact

assume you already know where the horizon is

never show how it emerges

Numerical relativity:

computes metrics

but does not expose global causal structure intuitively

Your project bridges this gap.

Vision statement:

A computational framework that takes spacetime data (analytic or numerical) and reconstructs its global causal structure — including event horizons — in real time, with rigorous mathematical grounding.

2. Core conceptual pillars (non-negotiable)

This tool is built around four deep GR concepts, not graphics:

Pillar I — Event horizons are global objects

Defined as 
∂
𝐽
−
(
𝐼
+
)
∂J
−
(I
+
)

Cannot be found locally

Must be computed retrospectively

Pillar II — Causality, not geometry, is fundamental

Distances are secondary

Light cones define structure

Horizons are causal boundaries

Pillar III — Horizons ≠ apparent horizons

Event horizons depend on the future

Apparent horizons depend on slicing

Your tool must show the difference

Pillar IV — Visualization must respect invariance

No misleading “space at a time” plots

Everything must be causal-structure faithful

This immediately elevates the project above standard “GR visualization”.

3. System architecture (high-level)

Think in three layers:

┌────────────────────────────────────┐
│  Visualization & Interaction Layer │
├────────────────────────────────────┤
│  Causal Structure Computation Core │
├────────────────────────────────────┤
│  Spacetime Input & Evolution Layer │
└────────────────────────────────────┘


Each layer is mathematically motivated.

4. Spacetime Input & Evolution Layer
What this layer accepts
(A) Analytic metrics

Schwarzschild

Vaidya (radiating collapse!)

Oppenheimer–Snyder

Kerr (later)

(B) Numerical spacetime data

Metric 
𝑔
𝜇
𝜈
(
𝑡
,
𝑥
)
g
μν
	​

(t,x)

From numerical relativity codes

Or synthetic collapse models

This layer defines:

(
𝑀
,
𝑔
𝜇
𝜈
)
(M,g
μν
	​

)

as raw spacetime data.

5. Causal Structure Computation Core (the heart)

This is where the project becomes research-grade.

5.1 Null geodesic engine

Numerically integrate:

𝑑
2
𝑥
𝜇
𝑑
𝜆
2
+
Γ
𝜈
𝜎
𝜇
𝑑
𝑥
𝜈
𝑑
𝜆
𝑑
𝑥
𝜎
𝑑
𝜆
=
0
(null)
dλ
2
d
2
x
μ
	​

+Γ
νσ
μ
	​

dλ
dx
ν
	​

dλ
dx
σ
	​

=0(null)

For:

large ensembles of null rays

launched from many spacetime points

in many directions

This gives you a causal flow field.

5.2 Numerical construction of 
𝐽
−
(
𝐼
+
)
J
−
(I
+
)

This is the key algorithmic insight.

Instead of guessing horizons:

Start with a grid of spacetime points 
𝑝
p

For each 
𝑝
p, launch null rays

Check whether any reach large 
𝑟
r at late times

Classify:

𝑝
∈
𝐽
−
(
𝐼
+
)
or not
p∈J
−
(I
+
)or not

Then:

𝐻
+
=
∂
𝐽
−
(
𝐼
+
)
H
+
=∂J
−
(I
+
)

This is definition-level faithful.

5.3 Retrodictive horizon construction

You then:

Identify null rays that asymptotically hover

Trace them backwards

Build the horizon surface

This makes the “retroactive” nature computationally explicit.

6. Visualization philosophy (this is crucial)
Rule 1: No misleading coordinates

You do not show raw Schwarzschild 
𝑡
,
𝑟
t,r unless explicitly requested.

Default views:

Penrose compactified

Conformal diagrams

Causal diamonds

Rule 2: Everything is interactive and falsifiable

User can:

Pick an event 
𝑝
p

Launch light rays

See whether it escapes or not

Change future evolution and watch horizon move

This kills confusion instantly.

Rule 3: Multiple simultaneous representations

Show the same causal structure in:

Schwarzschild

Kruskal

Penrose

Side-by-side.

This is pedagogically and scientifically powerful.

7. Research-level extensions (where papers live)
Extension A — Horizon vs apparent horizon mismatch

Let user:

Change spacetime slicing

Compute apparent horizons

Compare with true event horizon

This directly addresses:

“Horizon detection in numerical relativity”

Extension B — Dynamical horizons & trapping surfaces

Extend beyond event horizons:

Marginally trapped surfaces

Hayward horizons

Isolated horizons

This pushes into modern GR research.

Extension C — Quantum field theory input

Overlay:

mode propagation

near-horizon redshift

precursors to Hawking radiation

This ties directly to your QFT background.

Extension D — Causal set approximation

Discretize spacetime:

Nodes = events

Edges = causal relations

Compare continuum GR vs causal set theory.

This is very deep and modern.

8. Why this is genuinely new

Existing tools:

draw Penrose diagrams manually

visualize numerical data locally

do not compute global causal structure

Your tool:
✔ computes horizons from first principles
✔ exposes global causality dynamically
✔ makes nonlocality visible
✔ bridges GR, numerics, and QFT

This is publishable.

9. Execution roadmap (realistic)
Phase 1 — Prototype (1–2 months)

Schwarzschild + Vaidya

Null geodesic tracer

Penrose visualization

Phase 2 — Collapse physics (2–3 months)

Oppenheimer–Snyder

Horizon formation

Retroactive horizon visualization

Phase 3 — Research extensions (3–6 months)

Apparent vs event horizons

Dynamical horizons

Mode propagation

Phase 4 — Publication & release

Computational physics paper

Open-source platform

Interactive web demo

10. If I were your supervisor, I’d say this

“This project shows deep conceptual understanding, technical skill, and produces something genuinely useful to the GR community.”

That’s the level we’re talking about.