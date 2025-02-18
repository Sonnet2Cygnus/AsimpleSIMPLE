# AsimpleSIMPLE

This code implements a basic version of the 2D SIMPLE algorithm for solving simple incompressible fuid flow problems. 

## Features
- Modular structure with clear separation of functionalities
- Handles structured meshes for solving fluid flow problems
- Highly scalable

## Next Steps
- Expand to a 3D solver (partially completed)
- Develop a coupled flow-heat solver (compressible)
- Implement post-processing capabilities (.vtk, .dat)
- Extend to additional algorithms such as SIMPLEC, PIMPLE, and PISO

## References and Additional Information
- This codebase is under continuous development.
- The equations used in the code can be found in the accompanying Equation.pdf file.
- The core structure of this code is inspired by the Python implementation from the CFD 0to1 course. While it has been simplified, optimized, and extended, it largely retains the same framework as the original course material [1].
- Another key reference for this work is The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM® and Matlab®[2].

[1] https://www.bilibili.com/video/BV1Bo4y1s7NZ/?spm_id_from=333.337.search-card.all.click&vd_source=0ed3f9d7c4fbfbf5b11b4dbbe3c1fc20

[2] Moukalled, F., Mangani, L., & Darwish, M. (2016). The finite volume method in computational fluid dynamics: An advanced introduction with OpenFOAM® and Matlab®. Springer.
