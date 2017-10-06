! This program solves "Shallow Water Equations" in 1D using
! Richtmyer two-step Lax–Wendroff method with reference to John bukardt's code.
! 3 types of initial conditions and 3 boundary conditions are tested.
! IC: Sine wave, Dam break, Gaussian functions
! BC: Periodic, Free, reflective bcs
! @ author: VL Srinivas (https://www.linkedin.com/in/vl-srinivas/)

program main
implicit none

    integer, parameter:: nx =20 , nt = 250
    integer:: iter,i, j
    real(8) :: dx, dt, x_length, T_length
    character:: file_h, file_t, file_uh, file_x
    character(len=100)::fileplace, pico
    real(8), parameter:: g = 9.81d0, pi=3.141592653589793D0
    ! hm: cell centre values of h; h: nodal values of h;
    ! h_array: h values over the grid along with time
    real(8):: h(nx), h_array(nx,nt+1), hm(nx-1), x(nx)
    ! uhm: cell centre values of uh; uh: nodal values of uh;
    ! uh_array: uh values over the grid along with time
    real(8):: t(nt+1), uh(nx), uh_array(nx, nt+1), uhm(nx-1)
    ! RK4 variables
    real(8), dimension(nx):: kh1,kh2,kh3,kh4,kuh1,kuh2,kuh3,kuh4
    real:: t1, t2, t3       ! Time variables

    call cpu_time(t1)

    ! File names with extensions
    file_h = 'h'
    file_t = 't'
    file_uh = 'uh'
    file_x = 'x'

!   write(*,*) "Enter the following input values:"
    x_length = 1.0d0
    T_length = 0.50d0

    ! Creating spatial and time grid
    call linspace(nx, 0.0d0, x_length, X)
    !print *, 'X'
    !print *, X

    call linspace(nt, 0.0d0, T_length, T)
    !print *, 't'
    !print *, t

    ! Creating step size for time and space
    dx = x_length/(nx-1)
    !print *, 'dx', dx

    dt = T_length/(nt-1)
    print *, 'dt', dt

    call initial(pi, nx, X, h, uh, X_length)
    ! ----------------------- MAIN LOOP --------------------
    do iter = 1, nt

		! Lax-Wendroff's method for finite difference hyperbolic pde
    	! Taking half step at cell centres
        do i = 1, nx-1

            hm(i) =  ( h(i) + h(i+1) )/2.0d0 &
                - dt/2.0d0 *( uh(i+1) - uh(i) )/dx
            !print *, 'hm-',i,hm(i)

            uhm(i) = ( uh(i) + uh(i+1) )/2.0d0 - dt/2.0d0 * &
    			( ( uh(i+1)**2 /h(i+1) + g/2.0d0*h(i+1)**2 ) - &
    				( uh(i)**2 /h(i) + g/2.0d0 * h(i)**2) )

        end do
		! Taking the full time step at grid points defining the parameters
		! at the next time step only for internal grid points and excluding the
		! extreme points

        do i = 2, nx-1

            h(i) = h(i) - dt/dx *( uhm(i) - uhm(i-1) )
            !print *, 'h-',i,h(i)

            uh(i) = uh(i) - dt/dx *( uhm(i)**2/hm(i) + &
                0.5d0*g*hm(i)**2 - uhm(i-1)**2/hm(i-1) - &
                    0.5d0*g*hm(i-1)**2 )
        end do

		! ------------ Runge-Kutta 4th order scheme ----------------
		! RK step: 1
!		call func(dx, g, nx, h,				  uh, 				 kh1, kuh1)
!		! RK step: 2
!		call func(dx, g, nx, h+ dt/2.0d0, uh+ kuh1*dt/2.0d0, kh2, kuh2)
!		! RK step: 3
!		call func(dx, g, nx, h+ dt/2.0d0, uh+ kuh2*dt/2.0d0, kh3, kuh3)
!		! RK step: 4
!		call func(dx, g, nx, h + dt, 	  uh + kuh3, 		 kh4, kuh4)

!		! RK final step - adding all weights
!		h(1:nx) = h(1:nx) + ( kh1(1:nx) + 2.0d0*kh2(1:nx) + &
!			2.0d0*kh3(1:nx) + kh4(1:nx) ) * dt/6.0d0

!		uh(1:nx) = uh(1:nx) + ( kuh1(1:nx) + 2.0d0*kuh2(1:nx) + &
!			2.0d0*kuh3(1:nx) + kuh4(1:nx) ) * dt/6.0d0

        !print *, 'h-', h

!--------! Applying the boundary conditions
		call apply_BC(nx, h, uh, 3)

		! Entering the data into matrices with time
		h_array(1:nx, iter+1) = h(1:nx)

		uh_array(1:nx, iter+1) = uh(1:nx)

		! Writing the data to the files : ( h & uh )
		if ( mod(iter,1) == 0 )  then!call write_vec(file_h,fileplace,nx,h,iter)

            !write(fileplace,*) '/home/srinivas/Desktop/vls/files/'
            write(pico,*) iter

            open(unit=03, file= "/home/srinivas/Desktop/vls/files/h_"// &
            	trim(adjustl(pico))//".dat", status='unknown')

            write(03,*) "X steps:", nx
            write(03,*) "T steps:", nt

            !if (nx > 0) then
                ! Writing the data
                do j = 1, nx
                    write(03,100) real(X(j), kind=4), real(h(j), kind=4)
                end do
            !end if
            100 format(2x, 10(e12.6,' '))

            close(unit=03)
        end if
		!call write_vec(file_uh, fileplace, nx, h, iter)

	end do
	! ----------------------- END MAIN LOOP --------------------

	! Writing the arrays to files

	!call write_vec('T', fileplace, nt+1, T, 0)

	!call write_mat('h_array', fileplace, nx, nt+1, h_array)

	!call write_mat('uh_array', fileplace, nx, nt+1, uh_array)
	!call write_mat(filename, fileplace, x, y, matrix)

	call cpu_time(t2)

	t3 = t2 - t1

	write(*,*) "---------------------------------------------"
	write(*,*) "Shallow Water Equations 1D"
	write(*,*) "Spatial grid:", nx
	write(*,*) "Time grid:", nt
	write(*,*) "T_length:", t_length
	write(*,*) "X_length:", x_length
	write(*,*) "Spatial step size:", dx
	write(*,*) "Time step size:", dt
	write(*,*) "Time taken:", t3, 'sec'
	write(*,*) "---------------------------------------------"

 end program
!*******************************************************************************
subroutine linspace(n, first, last, A)
implicit none

    integer:: i
    integer, intent(in):: n
    real(8), intent(in):: first, last
    real(8), intent(out):: A(n)

    if (n==1) then
        a(1) = (first + last)/2.0d0

    else
        do i = 1, n

            A(i) = ( real(n-i, kind=8)*first + real(i-1, kind=8)*last )/ &
                   real(n-1,kind=8)
        end do

    end if

return
end subroutine
!*******************************************************************************
subroutine initial(pi, nx, X, h, uh, X_length)
implicit none

    integer, intent(in):: nx
	real(8), intent(in):: X(nx), pi, X_length
	real(8), intent(out):: h(nx), uh(nx)
	real(8)::a, b, c

    ! Sine wave initial condition
	!h(1:nx) = 2.0d0 + sin(2*pi* X(1:nx))

	! Dam break initial condition
	!h = 1.0d0
	!h(1:nx/2) = 3.0d0

    ! Gaussian function
    a=3.0d0
    b=X_length/2.0d0
    c=0.2d0
    h(1:nx) = a*exp(-(X(1:nx)-b)**2 /(2*c**2))

	uh(1:nx) = 0.0d0

return
end subroutine
!*******************************************************************************
subroutine apply_BC(nx, h, uh, bc)
implicit none

	integer, intent(in):: bc, nx
	real(8):: h(nx), uh(nx)

	! Applying Periodic boundary conditions
	if (bc==1) then
		h(1) = h(nx-1)
		h(nx) = h(2)
		uh(1) = uh(nx-1)
		uh(nx) = uh(2)

	! Applying Free Boundary Condition
	else if (bc==2) then
		h(1) = h(2)
		h(nx) = h(nx-1)
		uh(1) = uh(2)
		uh(nx) = uh(nx-1)

	! Applying Reflective BC for 'uh' & Free BC for 'h'
	else if (bc==3) then
		h(1) = h(2)
		h(nx) = h(nx-1)
		uh(1) = -uh(2)
		uh(nx) = -uh(nx-1)

	end if

return
end subroutine
!*******************************************************************************
subroutine write_vec(filename, fileplace, n, vector, iter)
! Write a vector of length n to a file
implicit none


	integer:: n, j, iter
	character:: filename, fileplace
	character(len=100):: pico
	real(8):: vector(n)

    write(pico,*) iter
	open(unit=01, file = "/home/srinivas/Desktop/vls/files/T_"//&
		trim(adjustl(pico))//".dat", status='unknown')

	write(01,*) "size:", sizeof(vector)

	if (n > 0) then
		! Writing the data
		do j = 1, n
			write(01,*) real(vector(j), kind=4)
		end do
	end if
	!100 format(2x, 10(e12.6,' '))

	close(unit=01)

return
end subroutine
!*******************************************************************************
subroutine write_mat(filename, fileplace, x, y, matrix)
! Writes a matrix of size x*y to a file with name 'filename'
! For every iteration there will be a fresh data file.
implicit none

	integer:: i, j, x, y
	real(8):: matrix(x,y)
	character:: filename, fileplace

	open(unit=02, file=trim(adjustl(fileplace))//trim(filename)//'.dat', &
		status='unknown')

	write(02,100) "size:", x, y

	if (x > 0 .and. y > 0) then
		! Writing the data
		do j = 1, y
			do i = 1, x
				write(01,100) matrix(i,j)
			end do
		end do

	end if
	100 format(2x, 10(e12.6,' '))

	close(unit=02)

return
end subroutine
!*******************************************************************************
subroutine func(dx, g, nx, hin, uhin, hout, uhout)
! Getting the function for time-marching for Runge-Kutta method
implicit none

	integer:: i, nx
	real(8), intent(in):: dx, g, hin(nx), uhin(nx)
	real(8), intent(out):: hout(nx), uhout(nx)


	! Applying Central difference scheme for spatial discretization

	do i = 1, nx

        print*, 'hin:', i,'=',hin(i)
        print*, 'uhin:', i,'=',uhin(i)


		hout(i) = -( uhin(i+1) - uhin(i-1) )/(2.0d0*dx)
		print *, 'hout:', hout(i)

		uhout(i) = -0.5d0/dx *( (uhin(i+1)**2 /hin(i+1) + 0.5d0*g*hin(i+1)**2) &
			- uhin(i+1)**2 /hin(i+1) + 0.5d0*g*hin(i+1**2))

	end do

return
end subroutine
!*******************************************************************************



















