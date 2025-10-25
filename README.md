# acPulse_vcLBM
This C++ code implements the Lattice Boltzmann Method based on an athermal D3Q19 velocity set and served as a foundation in a study of spurious acoustic reflexions created by the interaction of a Gaussian acoustic pulse and vertex-centered grid refinement interface, located in the middle of a quasi-2D domain, cf. https://doi.org/10.3390/fluids10020031.

<img width="2659" height="1560" alt="grafik" src="https://github.com/user-attachments/assets/8cbdb0f1-5742-4e65-b7dc-8ed3d0614d15" />
<br/>
<br/>
Included are several established collision operators:
  + Classical Bhatnagar-Gross-Krook model, cf. https://doi.org/10.1103/PhysRev.94.511
  + Multiple Relaxation Time model (Gram-Schmidt basis and raw-moments basis) , cf. https://doi.org/10.1103/PhysRevE.100.033305
  + Recursive-regularization model, cf. https://doi.org/10.48550/arXiv.1505.06900
  + Hybrid recursive-regularization model, cf. https://doi.org/10.1080/14685248.2018.1540879

... two vertex-centered grid refinement schemes:
  + "Standard" vertex-centered scheme with grid overlap, cf. https://doi.org/10.1006/jcph.1998.6089
  + More recent Direct-Coupling approach (no grid overlap), cf. https://doi.org/10.1016/j.jcp.2021.110667

... restriction operations for the information exchange during the fine-to-coarse communication step at the refinement interface:
  + Simple arithmetic average of fine non-equilibria, cf. https://doi.org/10.1016/j.jcp.2012.03.015
  + Anisotropic filtering of fine non-equilibria, cf. https://doi.org/10.1016/j.jcp.2013.07.037

... non-reflecting characteristic Dirichlet boundary conditions based on the local one-dimensional inviscid equations:
  + Classical version, cf. https://doi.org/10.1103/PhysRevE.78.046707
  + Including transverse terms, cf. http://dx.doi.org/10.1016/j.jcp.2015.08.044

Case parameters and model types have to be selected by the user within BoundaryConditions.txt by commenting/uncommenting the desired feature. The absolute path to this file needs to be set once within main (look for inputFile = "/absolute-path-to-BC-file/BoundaryConditions.txt";).
