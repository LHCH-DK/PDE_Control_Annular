# PDE_Control_Annular
    # To setup the optimal control problem, run the initialization script "Setup_main.m"
    # In "Setup_main.m" you can alter: 
    # Nvec- Number of Legendre modes
    # NN - Number of Fourier modes
    # M=NN/2 Number of trig modes
    # epsi - Tikonov Regularization paramter
    # a b - Inner and outer Radii of the Annulus
    # ua,ub - Lower and upper point-wise bound constraints
    # zd=@(r,theta) - The Desired state
    # f=@(r,theta) - Optimal Source term
    # IterMax - Max Number of SSN iterations
    # tolSSN - Tolerance SSN solver
    # tolKSP - Tolerance KSP solver
 
    # To run the code, execute the script "main.m". 
    # The script will solve the PDE-constrained optimal control problem specified by "Setup_main.m". 
    # The script prints the convergence history of the Krylov Subspace solver 
      and plots the optimal state and control.
