!     --------------------------------------------
    double precision function hesitation(rho)
!     --------------------------------------------
    implicit none

    double precision :: rho

    hesitation = 25*rho**0.2d0 / (1.d0-rho)**0.1d0

    return
    END function

!     --------------------------------------------
    double precision function dhes(rho)
!     --------------------------------------------
    implicit none

    double precision :: rho

    dhes = 5.d0*(2.d0-rho)/(2.d0*rho**0.8*(1.d0-rho)**1.1)

    return
    END function
    


    subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
               wave,s,amdq,apdq,num_aux)
! =====================================================

! This is an HLL Riemann solver for the ARZ traffic model.

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

! On output, wave contains the waves as jumps in q,
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
    double precision :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: rhoi, rhoim, qi, qim, hi, him, dhi, dhim
    double precision :: Q_hat_1, Q_hat_2, sl, sh, sr, rhoih, qih, hih, dhih


!     # split the jump in f at each interface into waves

    do i = 2-mbc, mx+mbc
        rhoi  = ql(1,i)
        rhoim = qr(1,i-1)
        rhoih = 0.5d0*(rhoi+rhoim)
        qi  = ql(2,i)
        qim = qr(2,i-1)
        qih = 0.5d0*(qi+qim)
        hi  = hesitation(rhoi)
        him = hesitation(rhoim)
        hih = hesitation(rhoih)
        dhi  = dhes(rhoi)
        dhim = dhes(rhoim)
        dhih = dhes(rhoih)

        !Compute speeds
        sl     = qim/rhoim - him - rhoim*dhim
        sh     = qih/rhoih - hih - rhoih*dhih
        s(1,i) = min(sl,sh)

        sh     = qih/rhoih - hih
        sr     = qi/rhoi - hi
        s(2,i) = max(sh,sr)

        if (s(1,i).gt.0.d0) then ! This is problematic
            write(*,*) 'what?'
            Q_hat_1 = rhoim
            Q_hat_2 = qim
            s(1,i) = 0.d0
        elseif (s(2,i).lt.0.d0) then ! This never happens
            write(*,*) 'whoa!'
            Q_hat_1 = rhoi
            Q_hat_2 = qi
            s(2,i) = 0.d0
        else
            if (abs(s(1,i)-s(2,i)).lt.1.e-10) then
                write(*,*) 'tiny!'
                Q_hat_1 = 0.5d0*(rhoi + rhoim)
                Q_hat_2 = 0.5d0*(qi + qim)
            else
                Q_hat_1 = ((qi-rhoi*hi) - (qim-rhoim*him) &
                            - s(2,i)*rhoi + s(1,i)*rhoim)/(s(1,i)-s(2,i))
                Q_hat_2 = ((qi**2.d0/rhoi - hi*qi) - (qim**2.d0/rhoim - him*qim) &
                            - s(2,i)*qi + s(1,i)*qim)/(s(1,i)-s(2,i))
            endif
        endif

        !if (abs(s(1,i)-s(2,i)).lt.1.e-10) then
        !    write(*,*) i, abs(s(1,i)-s(2,i))
        !endif

        ! Compute the waves
        wave(1,1,i) = Q_hat_1 - rhoim
        wave(2,1,i) = Q_hat_2 - qim
    
        wave(1,2,i) = rhoi - Q_hat_1
        wave(2,2,i) = qi   - Q_hat_2
        !write(*,*) i, wave(2,2,i)

    END DO


    ! Check signs here?
    do 220 m=1,meqn
        do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + wave(m,mw,i)
                endif
            90 END DO
    220 END DO

    return
    end subroutine rp1
