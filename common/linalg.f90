module linalg_mod
    implicit none
    
    contains

    subroutine lu_solve(N, A, b, x)
        ! Solves a general [A]x=b on an nxn matrix
        ! This replaces A (in place) with its LU decomposition (permuted row-wise)
      
          implicit none
      
          integer,intent(in) :: N
          real,dimension(N,N),intent(inout) :: A
          real,dimension(N),intent(in) :: b
          real,dimension(:),allocatable,intent(out) :: x
      
          integer,allocatable,dimension(:) :: indx
          integer :: D, info
      
          allocate(indx(N))
      
          ! Compute decomposition
          call lu_decomp(A, N, indx, D, info)
      
          ! if the matrix is nonsingular, then backsolve to find X
          if (info == 1) then
              write(*,*) 'Subroutine lu_decomp() failed. The given matrix is singular (i.e. no unique solution). Quitting...'
              stop
          else
              call lu_back_sub(A, N, indx, b, x)
          end if
      
          ! Cleanup
          deallocate(indx)
      
      end subroutine lu_solve
      
      
      !*******************************************************
      !*    LU decomposition routines used by test_lu.f90    *
      !*                                                     *
      !*                 F90 version by J-P Moreau, Paris    *
      !*    improved for F95 by Cory Goates, Logan, UT, USA  *
      !* --------------------------------------------------- *
      !* Reference:                                          *
      !*                                                     *
      !* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
      !*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
      !*  University Press, 1986" [BIBLI 08].                *
      !*                                                     *
      !*******************************************************
      
      
      subroutine lu_decomp(A, N, indx, D, code)
        ! Given an N x N matrix A, this routine replaces it by the LU
        ! decomposition of a rowwise permutation of itself. A and N  
        ! are input. indx is an output vector which records the row  
        ! permutation effected by the partial pivoting; D is output  
        ! as -1 or 1, depending on whether the number of row inter-  
        ! changes was even or odd, respectively. This routine is used
        ! in combination with lu_back_sub to solve linear equations or to 
        ! invert a matrix. Return code is 1 if matrix is singular.  
      
        implicit none
      
        real,dimension(N,N),intent(inout) :: A
        integer,intent(in) :: N
        integer,dimension(N),intent(out) :: indx
        integer,intent(out) :: code, D
      
        real,dimension(N) :: vv
        real,parameter :: tiny=1.5e-20
        integer :: i, j, k, imax
        real :: amax, dum, s
      
        ! Initialize
        D = 1
        code = 0
        imax = 0
      
        ! Loop over rows to get implicit scaling information
        do i=1,N
      
          ! Get largest element in this row
          amax=0.0
          do j=1,N
            if (abs(A(i,j)) > amax) then
              amax = abs(A(i,j))
            end if
          end do
      
          ! Check the largest element in this row is nonzero
          if (amax <= tiny) then
            code = 1 ! Singular matrix
            return
          end if
      
          ! Store scaling
          vv(i) = 1.0 / amax
      
        end do
      
        ! Loop over columns of Crout's method
        do j=1,N
      
          do i=1,j-1
            
            s = A(i,j)
      
            do k=1,i-1
              s = s - A(i,k)*A(k,j)
            end do
      
            A(i,j) = s
      
          end do
      
          ! Initialize search for largest pivot element
          amax = 0.0
          do i=j,N
        
            s = A(i,j)
            do k=1,j-1
              s = s - A(i,k)*A(k,j)
            end do
            A(i,j) = s
        
            ! Determine figure of merit for the pivot
            dum = vv(i)*abs(s)
            if (dum >= amax) then
              imax = i
              amax = dum
            end if
        
          end do
      
          ! Figure out if we need to interchange rows
          if (j /= imax) then
        
            ! Perform interchange
            do k=1,N
              dum = A(imax,k)
              A(imax,k) = A(j,k)
              A(j,k) = dum
            end do
        
            ! Update the sign of D since a row interchange has occurred
            D = -D
        
            ! Interchange the implicit scaling factor
            vv(imax) = vv(j)
        
          end if
      
          ! Store pivoting
          indx(j) = imax
      
          ! Divide by pivot element
          if (j /= N) then
            dum = 1.0 / A(j,j)
            do i=j+1,N
              A(i,j) = A(i,j)*dum
            end do
          end if
      
        end do
      
      end subroutine lu_decomp
      
      
      subroutine lu_back_sub(A, N, indx, b, x)
        ! Solves the set of N linear equations Ax = b.  Here A is     
        ! input, not as the matrix A but rather as its LU decomposition, 
        ! determined by the routine LUDCMP. indx is input as the permuta-
        ! tion vector returned by LUDCMP. b is input as the right-hand   
        ! side vector b. The solution vector is x. A, N, b and
        ! indx are not modified by this routine and can be used for suc- 
        ! cessive calls with different right-hand sides. This routine is 
        ! also efficient for plain matrix inversion.                     
      
        implicit none
      
        integer,intent(in) :: N
        real,dimension(N,N),intent(in) :: A
        real,dimension(N),intent(in) :: b
        integer,dimension(N),intent(in) :: indx
        real,dimension(:),allocatable,intent(out) :: x
      
        real :: sum
        integer :: ii,i,j,ll
      
        ! Initialize solution
        allocate(x, source=b)
      
        ! Set tracker to ignore leading zeros in b
        ii = 0
      
        ! Forward substitution
        do i=1,N
      
          ! Untangle pivoting
          ll = indx(i)
          sum = x(ll)
          x(ll) = x(i)
      
          ! If a nonzero element of b has already been encountered
          if (ii /= 0) then
            do J=ii,i-1
              sum = sum - A(i,J)*x(J)
            end do
      
          ! Check for first nonzero element of b
          else if(sum /= 0.0) then
            ii = i
          end if
      
          x(i) = sum
      
        end do
      
        ! Back substitution
        do i=N,1,-1
          sum = x(i)
          do j=i+1,N
            sum = sum - A(i,j)*x(j)
          end do
          x(i) = sum / A(i,i)
        end do
      
      end subroutine lu_back_sub
      


end module linalg_mod