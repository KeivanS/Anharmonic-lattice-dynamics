The developers should follow these conventions when adding their contribution to the project:
- **Documentation**: at the beginning of every fortran subroutine add a comment which explains what the subroutine does and why. The comments at the top of the subroutine should start with `!!`

    As an example:
    ```
    Subroutine svd(n,m,a,b,x)
    !! this subroutine performs a singular value decomposition to solve a set of overdetermined linear equations defined by the matrix a(n,m) and array b(n)
    !! it is used to solve for the force constants given force-displacement data
    ```
  
    In the above example, the first line explains what it does and the second line explains why it is used.

- **Style**: Standard fortran programming style should be adopted: 

    All subroutines must start with `IMPLICIT NONE`

    Intent of every variable must be specified (`IN`, `OUT`, or `INOUT`)

    `END` statements should be followed with the name of the corresponding subroutine
