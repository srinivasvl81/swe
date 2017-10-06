
! This program solves "Shallow Water Equations" in 2D using Richtmyer
! two-step Lax–Wendroff method with reference to John bukardt's 1D code.
! 3 types of initial conditions(IC) and 3 boundary conditions(BC) are tested.
! IC: Sine wave, Dam break, Gaussian functions
! BC: Periodic, Free, Reflective bcs
! @ author: VL Srinivas (https://www.linkedin.com/in/vl-srinivas/)

program main
implicit none

    integer, parameter:: nx=20, ny=20, nt=500
    real(8), parameter:: x_length=1.0d0, Y_length=1.0d0, T_length = 0.2d0
    integer:: iter, i, j
    real(8) :: dx, dy, dt, ck
    character(len=100):: fileplace, pico
    real(8), parameter:: g=9.81d0, pi=3.141592653589793D0
    real(8), dimension(nx,ny):: h, uh, vh
    real(8):: X(nx), Y(ny), t(nt+1)
    ! Arrays to store the values for each time step
    real(8), dimension(nx,ny,nt+1):: uh_array, vh_array, h_array
    ! Values at midpoints of grid points for half time step calculation
    real(8), dimension(nx-1,ny-1):: hmx, hmy, uhmx, uhmy, vhmx, vhmy
    ! RK4 variables
    real(8), dimension(nx):: kh1,kh2,kh3,kh4,kuh1,kuh2,kuh3,kuh4
    real:: t1, t2, t3       ! Time variables

    call cpu_time(t1)

!   write(*,*) "Enter the following input values:"
!    x_length = 1.0d0
!    T_length = 0.20d0

    ! Creating spatial and time grid
    call linspace(nx, 0.0d0, X_length, X)
    !print *, 'X'
    !print *, X
    call linspace(ny, 0.0d0, Y_length, Y)

    call linspace(nt, 0.0d0, T_length, T)
    !print *, 't'
    !print *, t

    ! Creating step size for time and space
    dx = X_length/(nx-1)
    !print *, 'dx', dx
    dy = Y_length/(ny-1)

    dt = T_length/(nt-1)

    call initial(pi, nx, ny, X, Y, h, uh, vh, X_length, Y_length)

    ck = max( maxval(h), abs(minval(h)) )

    ! ----------------------- MAIN LOOP --------------------
    do iter = 1, nt

        ! Solving using Ritchmyer 2-step lax-wendroff method
        call solver_lax_wendroff(nx,ny,dx,dy,dt,g,h,uh,vh)

!--------! Applying the boundary conditions
		call apply_BC(nx, ny, h, uh, vh, 3)

		! Entering the data into matrices with time
		h_array(1:nx, 1:ny, iter+1) = h(1:nx,1:ny)
		uh_array(1:nx, 1:ny, iter+1) = uh(1:nx,1:ny)
		vh_array(1:nx, 1:ny, iter+1) = vh(1:nx,1:ny)

		if ( iter == 1 )  then
            !call write_h(nx,ny,nt,dx,dy,dt,X,Y,h,iter)
        end if

		! Writing the data to the files : ( h & uh )
		if ( mod(iter,10) == 0 )  then
            call write_h(nx,ny,nt,dx,dy,dt,X,Y,h,iter)
        end if

		! Checking the time-step size
		!call check(dx,dt,nx,ny,h,uh,vh,ck)

        !if (iter==60) then


	end do
	! ----------------------- END MAIN LOOP --------------------

	! Writing the arrays to files

	call cpu_time(t2)

	t3 = t2 - t1

	write(*,*) "---------------------------------------------"
	write(*,*) "Shallow Water Equations 2D"
	write(*,*) "Spatial grid:", nx,'*',ny
	write(*,*) "Time grid:", nt
	write(*,*) "T_length:", real(t_length,kind=4)
	write(*,*) "X_length:", real(X_length,kind=4)
	write(*,*) "Y_length:", real(Y_length,kind=4)
	write(*,*) "Spatial step size-X:", real(dx,kind=4)
	write(*,*) "Spatial step size-Y:", real(dy,kind=4)
	write(*,*) "Time step size:", real(dt,kind=4)
	write(*,*) "Time taken:", t3, 'sec'
	write(*,*) "---------------------------------------------"
    print *, 'Extreme height',ck
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
subroutine initial(pi, nx, ny, X, Y, h, uh, vh, X_length, Y_length)
implicit none

    integer, intent(in):: nx, ny
	real(8), intent(in):: X(nx), Y(ny), pi, X_length, Y_length
	real(8), intent(out):: h(nx,ny), uh(nx,ny), vh(nx,ny)
	real(8)::a, x0, y0, sx, sy
	integer:: i, j

    ! Sine wave initial condition

	!h(1:nx, 1:ny) = 2.0d0 + sin(2*pi* X(1:nx))

	! Dam break initial condition
	!h = 1.0d0
	!h(1:nx/2, :) = 3.0d0

    ! Gaussian function
    a = 3.0d0
    x0 = X_length/2.0d0
    y0 = Y_length/2.0d0
    sx = 0.2d0
    sy = 0.2d0

    do j = 1, ny
    do i = 1, nx
        h(i,j) = a*exp(- ( (X(i)-x0)**2 /(2*sx**2) + (Y(j)-y0)**2 /(2*sy**2) ) )
    end do
    end do

	uh(1:nx, 1:ny) = 0.0d0
	vh(1:nx, 1:ny) = 0.0d0


return
end subroutine
!*******************************************************************************
subroutine apply_BC(nx, ny, h, uh, vh, bc)
implicit none

	integer, intent(in):: bc, nx, ny
	real(8), dimension(nx,ny):: h, uh, vh

	! Applying Periodic boundary conditions
	if (bc==1) then

		! South boundary
		h(:,1) = h(:,ny-1)
		uh(:,1) = uh(:,ny-1)
		vh(:,1) = vh(:,ny-1)

		! North boundary
		h(:,ny) = h(:,2)
		uh(:,ny) = uh(:,2)
		vh(:,ny) = vh(:,2)

		! West boundary
		h(1,:) = h(nx-1,:)
		uh(1,:) = uh(nx-1,:)
		vh(1,:) = vh(nx-1,:)

		! East boundary
		h(nx,:) = h(2,:)
		uh(nx,:) = uh(2,:)
		vh(nx,:) = vh(2,:)

	! Applying Free Boundary Condition
	else if (bc==2) then
		! South boundary
		h(:,1) = h(:,2)
		uh(:,1) = uh(:,2)
		vh(:,1) = vh(:,2)

		! North boundary
		h(:,ny) = h(:,ny-1)
		uh(:,ny) = uh(:,ny-1)
		vh(:,ny) = vh(:,ny-1)

		! West boundary
		h(1,:) = h(2,:)
		uh(1,:) = uh(2,:)
		vh(1,:) = vh(2,:)

		! East boundary
		h(nx,:) = h(nx-1,:)
		uh(nx,:) = uh(nx-1,:)
		vh(nx,:) = vh(nx-1,:)

	! Applying Reflective BC for 'uh' & Free BC for 'h'
	else if (bc==3) then
		! South boundary
		h(:,1) = h(:,2)
		uh(:,1) = uh(:,2)
		vh(:,1) = -vh(:,2)

		! North boundary
		h(:,ny) = h(:,ny-1)
		uh(:,ny) = uh(:,ny-1)
		vh(:,ny) = -vh(:,ny-1)

		! West boundary
		h(1,:) = h(2,:)
		uh(1,:) = -uh(2,:)
		vh(1,:) = vh(2,:)

		! East boundary
		h(nx,:) = h(nx-1,:)
		uh(nx,:) = -uh(nx-1,:)
		vh(nx,:) = vh(nx-1,:)

	end if

return
end subroutine
!*******************************************************************************
subroutine solver_lax_wendroff(nx,ny,dx,dy,dt,g,h,uh,vh)
implicit none

    integer:: i,j
	integer, intent(in):: nx, ny
	real(8), intent(in):: dx,dy,dt,g
	real(8), dimension(nx,ny):: h, uh, vh
	real(8), dimension(nx-1,ny-1):: hmx, hmy, uhmx, uhmy, vhmx, vhmy

	! Lax-Wendroff's method in 2D for finite difference hyperbolic pde
    ! Taking half step at cell centres
    do j = 1, ny-1
	do i = 1, nx-1

		! Half steps at mid-points for the continuity equation
        hmx(i,j) =  ( h(i,j) + h(i+1,j) )/2.0d0 - dt/2.0d0 *( (uh(i+1,j)-uh(i,j)) )/dx
	    hmy(i,j) =  ( h(i,j) + h(i,j+1) )/2.0d0 - dt/2.0d0 *( (vh(i,j+1)-vh(i,j)) )/dy
	    !print *, 'hm-',i,hm(i)

		! Half steps at mid-points for X-Momentum equation
		uhmx(i,j) = ( uh(i,j) + uh(i+1,j) )/2.0d0 - dt/(2.0d0*dx) * &
				( ( uh(i+1,j)**2 /h(i+1,j) + g/2.0d0*h(i+1,j)**2 ) - &
					( uh(i,j)**2 /h(i,j) + g/2.0d0*h(i,j)**2) )

		uhmy(i,j) = ( uh(i,j) + uh(i,j+1) )/2.0d0 - dt/(2.0d0*dy) * &
				( uh(i,j+1)*vh(i,j+1) /h(i,j+1) - uh(i,j)*vh(i,j) /h(i,j) )

		! Half steps at mid-points for Y-Momemtum equation
		vhmx(i,j) = ( vh(i+1,j) + vh(i,j) )/2.0d0 - dt/(2.0d0*dy) * &
				( uh(i+1,j)*vh(i+1,j) /h(i+1,j) - uh(i,j)*vh(i,j) /h(i,j) )

		vhmy(i,j) = ( vh(i,j) + vh(i+1,j) )/2.0d0 - dt/(2.0d0*dy) * &
				( ( vh(i,j+1)**2 /h(i,j+1) + g/2.0d0*h(i,j+1)**2 ) - &
					( vh(i,j)**2 /h(i,j) + g/2.0d0*h(i,j)**2) )

	end do
    end do

	! Taking the full time step at grid points defining the parameters
	! at the next time step only for internal grid points and excluding the
	! extreme points
    do j = 2, ny-1
    do i = 2, nx-1

		! Continuity equation
        h(i,j) = h(i,j) - dt/dx*( uhmx(i,j) - uhmx(i-1,j) ) - dt/dy*( vhmy(i,j) - vhmy(i,j-1) )
            !print *, 'h-',i,h(i)

  		! X-Momemtum equation
		uh(i,j) = uh(i,j) - dt/dx*( (uhmx(i,j)**2/hmx(i,j) + 0.5d0*g*hmx(i,j)**2) - &
            	(uhmx(i-1,j)**2/hmx(i-1,j) + 0.5d0*g*hmx(i-1,j)**2) ) - &
            	dt/dy *( uhmy(i,j)*vhmy(i,j)/hmy(i,j) - uhmy(i,j-1)*vhmy(i,j-1)/hmy(i,j-1))

        ! Y-Momemtum equation
        vh(i,j) = vh(i,j) - dt/dx*(uhmx(i,j)*vhmx(i,j)/hmx(i,j) - uhmx(i-1,j)*vhmx(i-1,j)/hmx(i-1,j)) &
            	- dt/dy*( (vhmy(i,j)**2/hmy(i,j) + 0.5d0*g*hmy(i,j)**2) - &
            	(vhmy(i,j-1)**2/hmy(i,j-1) + 0.5d0*g*hmy(i,j-1)**2) )

    end do
    end do

return
end subroutine
!*******************************************************************************
subroutine write_h(nx,ny,nt,dx,dy,dt,X,Y,h,iter)
implicit none

    integer:: i, j
    integer, intent(in):: iter,nx,ny,nt
    character(len=100):: pico
    real(8):: dx,dy,dt,X(nx),Y(ny),h(nx,ny)

    write(pico,*) iter

    open(unit=03, file= "/home/srinivas/Desktop/vls_2/files/h_"// &
            trim(adjustl(pico))//".dat", status='unknown')

    write(03,*) "X steps:", nx
    write(03,*) "Y steps:", ny
    write(03,*) "T steps:", nt
    write(03,*) "Time:", real(iter*dt,kind=4)
    write(03,*) "X step-size:", dx
    write(03,*) "Y step-size:", dx
    write(03,*) "T step-size:", dt
    write(03,*) "       X   ", "     Y    ", "    Height  "

    !if (nx > 0) then
    ! Writing the data
    do j = 1, ny
    do i = 1, nx

        write(03,100) real(X(i),kind=4), real(Y(j),kind=4), real(h(i,j),kind=4)

    end do
    end do
    !end if
    100 format(2x, 10(e12.6,' '))

    close(unit=03)

return
end subroutine
!*******************************************************************************
subroutine check(dx,dt,nx,ny,h,uh,vh,ck)
! This checks the time step size for cfl condition and
! sets the divergent height to the maximum from initial
implicit none

    integer:: i,j,nx,ny
    real(8):: dx, dt, u, v, cfl, ck
    real(8),dimension(nx,ny):: h,uh,vh,u1,v1

    do j = 1, ny
        do i = 1, nx

            u1(i,j) = uh(i,j)/h(i,j)
            v1(i,j) = vh(i,j)/h(i,j)

            ! Setting the divergent height value to maximum
            if (h(i,j) > ck) h(i,j) = ck

            if (h(i,j) < -ck) h(i,j) = -ck

        end do
    end do

    u = max( maxval(u1), abs(minval(u1)) )
    v = max( maxval(v1), abs(minval(v1)) )

    cfl = dx/( sqrt(u**2 + v**2)*sqrt(2.0) )

    if (dt > cfl) dt = cfl/2.0d0

return
end subroutine
!*******************************************************************************











