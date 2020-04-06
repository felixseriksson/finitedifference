# extend from 1d to 2d: apply multivariable calculus rules, build 2d grid and discretize
# by keeping x, then y constant
# uij = u(xi, yj)
# ui+1j+1 = uij + (dx*d/dx + dy*d/dy)(uij) + hot (simple taylor series expansion)
# can do FD, BD and CD in x, y and time
# du/dt + cdu/dx + cdu/dy = 0
# FD in time, BD in x and y:
# uij^n+1 = uij^n - cdt/dx(uij^n - ui-1j^n) - cdt/dy(uij^n-uij-1^n)
# ICs/BCs:
# square wave but in 3d, so cube (?) wave? hypersquare wave ?