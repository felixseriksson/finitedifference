# FD in time, BD in space (as with 2d burgers), CD for pressure (usual difference for laplacian)
# discretized:
# LHS = (uij^n+1 - uij^n)/dt + uij^n((uij^n - ui-1j^n)/dx) + vij^n((uij^n - uij-1^n)/dy) = 
# RHS = -1/rho((pi+1j^n - pi-1j^n)/(2dx)) + nu((ui+1j^n - 2uij^n + ui-1j^n)/dx2 + (uij+1^n - 2uij^n + uij-1^n)/dy2)
# same for v also...
# poisson:
# LHS = (pi+1j^n - pij^n + pi-1j^n)/dx2 + (pij+1^n - pij^n + pij-1^n)/dy2 = 
# RHS = -rho(((ui+1j^n - ui-1j^n)/(2dx))^2 + 2((uij+1^n - uij-1^n)/(2dy))((vi+1j^n - vi-1j^n)/(2dx)) + ((vij+1^n - vij-1^n)/(2dy))^2) + rho/dt((ui+1j^n - ui-1j^n)/(2dx) + (vij+1^n - vij-1^n)/(2dy))
# transpose by isolating, uij^n+1, vij^n+1 and pij^n
# ICs:
# u, v, p = 0 everywehre
# BCs:
# u = 1 at y = 2
# u, v = 0 at x = 0, 2 and y = 0 (basically no-slip everywhere except at y = 2)
# p = 0 at y = 2
# dp/dy = 0 at y = 0 and dp/dx = 0 at x = 0, 2
# we will have actual time steps - nt - and artificial time steps - nit - for iteration of pressure inside outer loop