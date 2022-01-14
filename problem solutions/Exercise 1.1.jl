### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ ba3695dc-74d2-11ec-1c01-453e48962a5d
md"""
## Exercise 1.1
"""

# ╔═╡ 61046e8c-239f-4173-bade-9fe5bcb3b0e1
md"""
### 1. 
```math
\begin{align}
\int_{-\infty}^{\infty}\mathrm{d}x~e^{-ax^2}
= & \sqrt{\left(\int_{-\infty}^{\infty}\mathrm{d}x~e^{-ax^2}\right)^2}
=  \sqrt{\left(\int_{-\infty}^{\infty}\mathrm{d}x~e^{-ax^2}\right)\left(\int_{-\infty}^{\infty}\mathrm{d}y~e^{-ay^2}\right)}\\
= & \sqrt{\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}\mathrm{d}x\mathrm{d}y~e^{-a(x^2+y^2)}}
= \sqrt{\int_{0}^{\infty}\mathrm{d}r\int_{0}^{2\pi}\mathrm{d}\phi~re^{-ar^2}}\\
= & \sqrt{2\pi\int_{0}^{\infty}\mathrm{d}r~re^{-ar^2}}
= \sqrt{\pi\int_{0}^{\infty}\mathrm{d}r^2~e^{-ar^2}}\\
= & \sqrt{-\frac{\pi}{a}e^{-ar^2}\bigg|_0^\infty}
= \sqrt{\frac{\pi}{a}}\\
\int_{-\infty}^{\infty}\mathrm{d}x e^{-(x-\mu)^2/(2\sigma^2)}
= & \int_{-\infty}^{\infty}\mathrm{d}x~e^{-x^2/(2\sigma^2)}=\sqrt{2\pi\sigma^2}
\end{align}
```

### 2. 
```math
\begin{align}
\int_{-\infty}^{\infty}\mathrm{d}x~x^2e^{-ax^2} &= -\int_{-\infty}^{\infty}\mathrm{d}x~\frac{\mathrm{d}}{\mathrm{d}a}e^{-ax^2} = -\frac{\mathrm{d}}{\mathrm{d}a}\int_{-\infty}^{\infty}\mathrm{d}x~e^{-ax^2}= -\frac{\mathrm{d}}{\mathrm{d}a}\sqrt{\frac{\pi}{a}}= \frac{1}{2}\sqrt{\frac{\pi}{a^3}}\\
\int_{-\infty}^{\infty}\mathrm{d}x~xe^{-ax^2} &= \frac{1}{2}\int_{-\infty}^{\infty}\mathrm{d}x^2~e^{-ax^2} = \frac{1}{2a}
\end{align}
```

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─ba3695dc-74d2-11ec-1c01-453e48962a5d
# ╠═61046e8c-239f-4173-bade-9fe5bcb3b0e1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
