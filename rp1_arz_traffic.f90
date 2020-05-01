!     --------------------------------------------
double precision function hesitation(rho)
!     --------------------------------------------
    implicit none

    double precision :: rho

    hesitation = 25*rho**0.2d0 / (1.d0-rho)**0.1d0

    return
end function

!     --------------------------------------------
double precision function dhes(rho)
!     --------------------------------------------
    implicit none

    double precision :: rho

    dhes = 5.d0*(2.d0-rho)/(2.d0*rho**0.8*(1.d0-rho)**1.1)

    return
end function
    

! =====================================================
subroutine rp1(maxm,meqn,mwaves,num_aux,mbc,mx,ql,qr,auxl,auxr, &
               fwave,s,amdq,apdq)
! =====================================================

! This is an f-wave Riemann solver for the ARZ traffic model.

! Conservative form of the equations:

!       \rho_t + (q-\rho h(\rho))_x    = 0
!       q_t + (q^2/\rho - q h(\rho))_x = (\rho (U(\rho)+h(\rho)) - q)/tau

! Variable meanings:
!       \rho : density
!       u    : velocity
!       q = \rho (u + h(\rho))
!       U    : desired velocity
!       h    : hesitation
!       \tau : relaxation time
!
! On input, ql contains the state vector at the left edge of each cell
! qr contains the state vector at the right edge of each cell

! On output, fwave contains the waves as jumps in f,
! s the speeds,
! 
! amdq = A^- Delta q,
! apdq = A^+ Delta q,
! the decomposition of the flux difference
! f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.

! Note that the ith Riemann problem has left state qr(:,i-1)
!                                  and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit none

    double precision :: hesitation, dhes
    integer :: num_aux, mbc, maxm, meqn, mwaves, mx, m, mw
    integer :: i
    double precision :: auxl(num_aux,1-mbc:maxm+mbc)
    double precision :: auxr(num_aux,1-mbc:maxm+mbc)
    double precision :: fwave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: rhoi, rhoim, qi, qim, hi, him, dhi, dhim, denom
    double precision :: df1, df2, b1, b2

    do i = 2-mbc, mx+mbc
        rhoi  = ql(1,i)
        rhoim = qr(1,i-1)
        qi  = ql(2,i)
        qim = qr(2,i-1)
        hi  = hesitation(rhoi)
        him = hesitation(rhoim)
        dhi  = dhes(rhoi)
        dhim = dhes(rhoim)

        ! Compute jump in flux
        df1 = (qi-rhoi*hi) - (qim-rhoim*him)
        df2 = (qi**2.d0/rhoi - hi*qi) - (qim**2.d0/rhoim - him*qim)

        s(1,i) = qim/rhoim - him - rhoim*dhim
        s(2,i) = qi/rhoi - hi

        if (s(1,i).le.0.d0) then ! 1 wave going each direction
            ! Decompose flux jump in terms of flux jacobian eigenvectors
            denom = qi/rhoi - qim/rhoim + rhoi*dhi
            b1 = ( (qi/rhoi + rhoi*dhi)*df1 - df2 ) / denom
            b2 = ( -qim/rhoim*df1 + df2 ) / denom

            ! Compute the waves.
            fwave(1,1,i) = b1
            fwave(2,1,i) = b1 * qim/rhoim
            fwave(1,2,i) = b2
            fwave(2,2,i) = b2*( qi/rhoi + rhoi*dhi)
        else ! 2 right-going waves
            s(1,i) = qi/rhoi - hi - rhoi*dhi
            denom = rhoi*dhi
            b1 = ( (qi/rhoi + rhoi*dhi)*df1 - df2 ) / denom
            b2 = ( -qi/rhoi*df1 + df2 ) / denom
            fwave(1,1,i) = b1
            fwave(2,1,i) = b1 * qi/rhoi
            fwave(1,2,i) = b2
            fwave(2,2,i) = b2*( qi/rhoi + rhoi*dhi)
            if (s(1,i).le.0.d0) then
                ! Transonic 1-Shock
                s(1,i) = (qim-qi + rhoi*hi - rhoim*him)/(rhoim-rhoi)
            endif
        endif
    end do


    ! compute the leftgoing and rightgoing fluctuations:
    ! Note that here we assume s(1,i) < 0   and   s(2,i) > 0.
    do m=1,meqn
        do i = 2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + fwave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + fwave(m,mw,i)
                endif
            end do
        end do
    end do

    return
    end subroutine rp1



