### Continuum Puyo

This code is trying to solve the PDE

$$ \frac{\partial h}{\partial t} = \left(1-\frac2N\right) + \frac{a^2p^2}2\nabla^2h - \frac{c}N\left(\nabla h\right)^2 - ap^2\lvert\nabla h \rvert + \eta(x, t) $$

The main difference from KPZ is the $\lvert\nabla h\rvert$ term, so I also play around with giving it a strength $\lambda_\text{abs}$.

