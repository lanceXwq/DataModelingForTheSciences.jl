### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ adf97b68-9c03-11ec-2693-a95698a75d57
md"
We start from rewriting $\alpha$:
```math
\begin{align}
\alpha_r(r')
=\min\left(1,A_r(r')\right)=\frac{1}{2}\left[1+A_r(r')-\left|1-A_r(r')\right|\right].
\end{align}
```
Then
```math
\begin{align}
\pi(r)\alpha_r(r')Q_r(r')
&=\frac{\pi(r)Q_r(r')}{2}\left[1+A_r(r')-\left|1-A_r(r')\right|\right]\\
&=\frac{\pi(r)Q_r(r')}{2}\left[1+\frac{\pi(r')Q_{r'}(r)}{\pi(r)Q_r(r')}-\left|1-\frac{\pi(r')Q_{r'}(r)}{\pi(r)Q_r(r')}\right|\right]\\
&=\frac{1}{2}\left[\pi(r)Q_r(r')+\pi(r')Q_{r'}(r)-\left|\pi(r)Q_r(r')-\pi(r')Q_{r'}(r)\right|\right]\\
&=\min\left(\pi(r)Q_r(r'),\pi(r')Q_{r'}(r)\right).
\end{align}
```
Swapping $r$ and $r'$ yields $\pi(r')\alpha_{r'}(r)Q_{r'}(r)=\min\left(\pi(r')Q_{r'}(r),\pi(r)Q_r(r')\right)$, which is clearly equivalent to $\pi(r')\alpha_{r'}(r)Q_{r'}(r)$.
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─adf97b68-9c03-11ec-2693-a95698a75d57
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
