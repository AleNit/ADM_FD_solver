# ADM_FD_solver
Matlab scripts for solving the steady, variable-coefficient Advection-Diffusion-Migration (ADM) equation by a finite difference approximation in 1D, 2D or 3D in a cylindrical coordinate system over non-uniform structured grids. The discretized equation takes the form:

$$ A \nabla^2 q + B \nabla q - C q = r \, $$ 

where $$A$$, $$B$$ and $$C$$ are non-constant scalar fields.
Derivatives are approximated by a 2-nd order accurate centered difference over a cell-centered grid. Grid stretching can be prescribed. Dirichelet/Neumann/periodic boundary conditions are prescribed by a ghost cell method. The generic boundary condition is formulated as a Robin condition. The linear system arising from 1D and 2D examples is solved by a direct method, whereas that arising from the 3D example is solved by the BICGSTAB algorithm.

The discretization arrangement follows the research papers:

    Farnell, L. "Solution of Poisson equations on a nonuniform grid." Journal of Computational Physics 35.3 (1980): 408-425.
    
