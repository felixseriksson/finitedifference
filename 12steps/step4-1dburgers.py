# 1D Burgers' equation
# du/dt + u*du/dx = mu*d2u/dx2
# FD in time, BD in first order space, CD in 2nd order space
# ui^n+1 = ui^n - ui^n(dt/dx)(ui^n - ui-1^n) + nu(dt/dx2)(ui+1^n - 2*ui^n + ui-1^n)