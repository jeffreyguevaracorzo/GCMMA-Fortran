module GCMMA_Routines
    implicit none
    ! Base routines
    interface Zeros
        module procedure ZerosMatrixDP
        module procedure ZerosVectorDP
    end interface Zeros
    interface Ones
        module procedure OnesMatrixDP
        module procedure OnesVectorDP
    end interface Ones
    interface MultMatDiag
        module procedure MultMatDiagReal
        module procedure MultMatDiagDP
    end interface MultMatDiag
    interface Diag
        module procedure DiagMatrixReal
        module procedure DiagMatrixDP
    end interface Diag
    ! additional interfaces for [A]{x}={B}
    interface Solver
        module procedure SolverReal
        module procedure SolverDP
    end interface Solver
    ! GCMMA routines
    interface GCMMA_Sub
        module procedure GCMMAsubReal
        module procedure GCMMAsubDP
    end interface GCMMA_Sub
    interface Subsolve
        module procedure SubSolveReal
        module procedure SubSolveDP
    end interface Subsolve
    interface KKTCheck
        module procedure KKTCheckReal
        module procedure KKTCheckDP
    end interface KKTCheck
    interface Concheck
        module procedure ConcheckReal
        module procedure ConcheckDP
    end interface Concheck
    interface raaUpdate
        module procedure raaUpdateReal
        module procedure raaUpdateDP
    end interface raaUpdate
    interface Asymp
        module procedure AsympReal
        module procedure AsympDP
    end interface Asymp
contains
    ! zeros
    function ZerosMatrixDP(row,col) Result(Array)
        implicit none
        integer, intent(in)                             :: row
        integer, intent(in)                             :: col
        double precision, dimension(:,:), allocatable   :: Array
        allocate(Array(row,col))
        Array = 0.0d0
    end function ZerosMatrixDP
    
    function ZerosVectorDP(row) Result(Array)
        implicit none
        integer, intent(in)                             :: row
        double precision, dimension(:), allocatable     :: Vector
        double precision, dimension(:), allocatable     :: Array
        allocate(Vector(row))
        Vector = 0.0d0
        Array = Vector
    end function ZerosVectorDP
    
    ! ones
    function OnesMatrixDP(row,col) Result(Array)
        implicit none
        integer, intent(in)                             :: row
        integer, intent(in)                             :: col
        double precision, dimension(:,:), allocatable   :: Array
        allocate(Array(row,col))
        Array = 1.0d0
    end function OnesMatrixDP

    function OnesVectorDP(row) Result(Array)
        implicit none
        integer, intent(in)                             :: row
        double precision, dimension(:), allocatable     :: Array
        allocate(Array(row))
        Array = 1.0d0
    end function OnesVectorDP

    ! MultMatDiag
    function MultMatDiagReal(ArrayIn,VectorIn) Result(ArrayOut)
        implicit none
        integer                                                 :: i
        real, dimension(:,:), allocatable, intent(in)           :: VectorIn
        real, dimension(:,:), allocatable, intent(in)           :: ArrayIn
        real, dimension(:,:), allocatable                       :: ArrayOut
        allocate(ArrayOut(size(ArrayIn,1),size(ArrayIn,2)))
        do i = 1, size(ArrayIn,2), 1
            ArrayOut(:,i) = VectorIn(i,1)*ArrayIn(:,i)
        end do
    end function MultMatDiagReal

    function MultMatDiagDP(ArrayIn,VectorIn) Result(ArrayOut)
        implicit none
        integer                                                     :: i
        double precision, dimension(:,:), allocatable, intent(in)   :: VectorIn
        double precision, dimension(:,:), allocatable, intent(in)   :: ArrayIn
        double precision, dimension(:,:), allocatable               :: ArrayOut
        allocate(ArrayOut(size(ArrayIn,1),size(ArrayIn,2)))
        do i = 1, size(ArrayIn,2), 1
            ArrayOut(:,i) = VectorIn(i,1)*ArrayIn(:,i)
        end do
    end function MultMatDiagDP

    ! Diagonal Matrix
    function DiagMatrixReal(ArrayIn) Result(ArrayOut)
        implicit none
        integer                                         :: i
        real, dimension(:,:), allocatable, intent(in)   :: ArrayIn
        real, dimension(:,:), allocatable               :: ArrayOut
        i = maxval([size(ArrayIn,1),size(ArrayIn,2)])
        allocate(ArrayOut(i,i))
        ArrayOut = 0.0d0
        if (size(ArrayIn,1).gt.1) then
            do i = 1, size(ArrayIn,1), 1
                ArrayOut(i,i) = ArrayIn(i,1)
            end do
        else
            do i = 1, size(ArrayIn,2), 1
                ArrayOut(i,i) = ArrayIn(1,i)
            end do
        end if 
    end function DiagMatrixReal

    function DiagMatrixDP(ArrayIn) Result(ArrayOut)
        implicit none
        integer                                                     :: i
        double precision, dimension(:,:), allocatable, intent(in)   :: ArrayIn
        double precision, dimension(:,:), allocatable               :: ArrayOut
        i = maxval([size(ArrayIn,1),size(ArrayIn,2)])
        allocate(ArrayOut(i,i))
        ArrayOut = 0.0d0
        if (size(ArrayIn,1).gt.1) then
            do i = 1, size(ArrayIn,1), 1
                ArrayOut(i,i) = ArrayIn(i,1)
            end do
        else
            do i = 1, size(ArrayIn,2), 1
                ArrayOut(i,i) = ArrayIn(1,i)
            end do
        end if 
    end function DiagMatrixDP

    ! GCMMA Routine
    subroutine GCMMAsubReal(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp, &                                
                        m,n,iter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d)      
        implicit none
        ! in/OutputArguments
        real, dimension(:,:), allocatable, intent(out)  :: xmma
        real, dimension(:,:), allocatable, intent(out)  :: ymma
        real, intent(out)                               :: zmma
        real, dimension(:,:), allocatable, intent(out)  :: lam
        real, dimension(:,:), allocatable, intent(out)  :: xsi
        real, dimension(:,:), allocatable, intent(out)  :: eta
        real, dimension(:,:), allocatable, intent(out)  :: mu
        real, intent(out)                               :: zet
        real, dimension(:,:), allocatable, intent(out)  :: s
        real, intent(out)                               :: f0app
        real, dimension(:,:), allocatable, intent(out)  :: fapp
        integer, intent(in)                             :: m
        integer, intent(in)                             :: n
        integer, intent(in)                             :: iter
        real, intent(in)                                :: epsimin
        real, intent(in)                                :: a0
        real, intent(in)                                :: raa0
        real, intent(in)                                :: f0val
        real, dimension(:,:), allocatable, intent(in)   :: fval
        real, dimension(:,:), allocatable, intent(in)   :: xval
        real, dimension(:,:), allocatable, intent(in)   :: xmin
        real, dimension(:,:), allocatable, intent(in)   :: xmax
        real, dimension(:,:), allocatable, intent(in)   :: raa
        real, dimension(:,:), allocatable, intent(in)   :: df0dx
        real, dimension(:,:), allocatable, intent(in)   :: dfdx
        real, dimension(:,:), allocatable, intent(inout):: low
        real, dimension(:,:), allocatable, intent(inout):: upp
        real, dimension(:,:), allocatable, intent(in)   :: a
        real, dimension(:,:), allocatable, intent(in)   :: c
        real, dimension(:,:), allocatable, intent(in)   :: d
        ! InternalArguments
        integer                                         :: i,j
        real                                            :: move
        real                                            :: albefa
        real                                            :: asyinit
        real                                            :: asyincr
        real                                            :: asydecr
        real, dimension(:,:), allocatable               :: eeen
        real, dimension(:,:), allocatable               :: eeem
        real, dimension(:,:), allocatable               :: zeron
        real, dimension(:,:), allocatable               :: zzz
        real, dimension(:,:), allocatable               :: zzz1
        real, dimension(:,:), allocatable               :: zzz2
        real, dimension(:,:), allocatable               :: factor
        real, dimension(:,:), allocatable               :: lowmin
        real, dimension(:,:), allocatable               :: lowmax
        real, dimension(:,:), allocatable               :: uppmin
        real, dimension(:,:), allocatable               :: uppmax
        real, dimension(:,:), allocatable               :: alfa
        real, dimension(:,:), allocatable               :: beta
        real, dimension(:,:), allocatable               :: xmami
        real, dimension(:,:), allocatable               :: xmamieps
        real, dimension(:,:), allocatable               :: xmamiinv
        real, dimension(:,:), allocatable               :: ux1
        real, dimension(:,:), allocatable               :: ux2
        real, dimension(:,:), allocatable               :: xl1
        real, dimension(:,:), allocatable               :: xl2
        real, dimension(:,:), allocatable               :: uxinv
        real, dimension(:,:), allocatable               :: xlinv
        real, dimension(:,:), allocatable               :: p0
        real, dimension(:,:), allocatable               :: q0
        real, dimension(:,:), allocatable               :: pq0
        real, dimension(:,:), allocatable               :: P
        real, dimension(:,:), allocatable               :: Q
        real, dimension(:,:), allocatable               :: PQ
        real, dimension(:,:), allocatable               :: b
        real, dimension(:,:), allocatable               :: r
        real                                            :: r0
        ! PROCEDURE
        eeen = ones(n,1)
        zeron = zeros(n,1)
        albefa = 0.1
        move = 0.5

        ! Calculation of the bounds alfa and beta
        zzz1 = low + albefa*(xval-low)
        zzz2 = xval - move*(xmax-xmin)
        zzz  = max(zzz1,zzz2)
        alfa = max(zzz,xmin)
        zzz1 = upp - albefa*(upp-xval)
        zzz2 = xval + move*(xmax-xmin)
        zzz  = min(zzz1,zzz2)
        beta = min(zzz,xmax)

        ! Calculations of p0, q0, P, Q, r and b
        xmami = xmax - xmin
        xmamieps = 0.00001*eeen
        xmami = max(xmami,xmamieps)
        xmamiinv = eeen/xmami
        ux1 = upp - xval
        ux2 = ux1*ux1
        xl1 = xval - low
        xl2 = xl1*xl1
        uxinv = eeen/ux1
        xlinv = eeen/xl1
        p0 = zeron
        q0 = zeron
        p0 = max(df0dx,0.0)
        q0 = max(-1.0*df0dx,0.0)
        pq0 = p0 + q0
        p0 = p0 + 0.001*pq0
        q0 = q0 + 0.001*pq0
        p0 = p0 + raa0*xmamiinv
        q0 = q0 + raa0*xmamiinv
        p0 = p0*ux2
        q0 = q0*xl2
        r0 = f0val - maxval(matmul(transpose(p0),uxinv)) - maxval(matmul(transpose(q0),xlinv))
        !P = zeros(m,n)
        !Q = zeros(m,n)
        P = max(dfdx,0.0)
        Q = max(-1.0*dfdx,0.0)
        PQ = P + Q
        P = P + 0.001*PQ
        Q = Q + 0.001*PQ
        P = P + matmul(raa,transpose(xmamiinv))
        Q = Q + matmul(raa,transpose(xmamiinv))
        P = MultMatDiag(P,ux2)
        Q = MultMatDiag(Q,xl2)
        r = fval - matmul(P,uxinv) - matmul(Q,xlinv)
        b = -r

        ! Solving the subproblem by a primal-dual Newton method
        call Subsolve(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,m,n,epsimin,low,upp,alfa,beta, &
                    p0,q0,P,Q,a0,a,b,c,d)

        ! Calculations of f0app and fapp.
        ux1 = upp - xmma
        xl1 = xmma - low
        uxinv = eeen/ux1
        xlinv = eeen/xl1
        f0app = r0 + maxval(matmul(transpose(p0),uxinv)) + maxval(matmul(transpose(q0),xlinv))
        fapp  = r + matmul(P,uxinv) + matmul(Q,xlinv)
    end subroutine GCMMAsubReal

    subroutine GCMMAsubDP(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp, &                                
                        m,n,iter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d)    
        implicit none
        ! In/OutputArguments
        double precision, dimension(:,:), allocatable, intent(out)  :: xmma
        double precision, dimension(:,:), allocatable, intent(out)  :: ymma
        double precision, intent(out)                               :: zmma
        double precision, dimension(:,:), allocatable, intent(out)  :: lam
        double precision, dimension(:,:), allocatable, intent(out)  :: xsi
        double precision, dimension(:,:), allocatable, intent(out)  :: eta
        double precision, dimension(:,:), allocatable, intent(out)  :: mu
        double precision, intent(out)                               :: zet
        double precision, dimension(:,:), allocatable, intent(out)  :: s
        double precision, intent(out)                               :: f0app
        double precision, dimension(:,:), allocatable, intent(out)  :: fapp
        integer, intent(in)                                         :: m
        integer, intent(in)                                         :: n
        integer, intent(in)                                         :: iter
        double precision, intent(in)                                :: epsimin
        double precision, intent(in)                                :: a0
        double precision, intent(in)                                :: raa0
        double precision, intent(in)                                :: f0val
        double precision, dimension(:,:), allocatable, intent(in)   :: fval
        double precision, dimension(:,:), allocatable, intent(in)   :: xval
        double precision, dimension(:,:), allocatable, intent(in)   :: xmin
        double precision, dimension(:,:), allocatable, intent(in)   :: xmax
        double precision, dimension(:,:), allocatable, intent(in)   :: raa
        double precision, dimension(:,:), allocatable, intent(in)   :: df0dx
        double precision, dimension(:,:), allocatable, intent(in)   :: dfdx
        double precision, dimension(:,:), allocatable, intent(inout):: low
        double precision, dimension(:,:), allocatable, intent(inout):: upp
        double precision, dimension(:,:), allocatable, intent(in)   :: a
        double precision, dimension(:,:), allocatable, intent(in)   :: c
        double precision, dimension(:,:), allocatable, intent(in)   :: d
        ! InternalArguments
        integer                                                     :: i,j
        double precision                                            :: move
        double precision                                            :: albefa
        double precision                                            :: asyinit
        double precision                                            :: asyincr
        double precision                                            :: asydecr
        double precision, dimension(:,:), allocatable               :: eeen
        double precision, dimension(:,:), allocatable               :: eeem
        double precision, dimension(:,:), allocatable               :: zeron
        double precision, dimension(:,:), allocatable               :: zzz
        double precision, dimension(:,:), allocatable               :: zzz1
        double precision, dimension(:,:), allocatable               :: zzz2
        double precision, dimension(:,:), allocatable               :: factor
        double precision, dimension(:,:), allocatable               :: lowmin
        double precision, dimension(:,:), allocatable               :: lowmax
        double precision, dimension(:,:), allocatable               :: uppmin
        double precision, dimension(:,:), allocatable               :: uppmax
        double precision, dimension(:,:), allocatable               :: alfa
        double precision, dimension(:,:), allocatable               :: beta
        double precision, dimension(:,:), allocatable               :: xmami
        double precision, dimension(:,:), allocatable               :: xmamieps
        double precision, dimension(:,:), allocatable               :: xmamiinv
        double precision, dimension(:,:), allocatable               :: ux1
        double precision, dimension(:,:), allocatable               :: ux2
        double precision, dimension(:,:), allocatable               :: xl1
        double precision, dimension(:,:), allocatable               :: xl2
        double precision, dimension(:,:), allocatable               :: uxinv
        double precision, dimension(:,:), allocatable               :: xlinv
        double precision, dimension(:,:), allocatable               :: p0
        double precision, dimension(:,:), allocatable               :: q0
        double precision, dimension(:,:), allocatable               :: pq0
        double precision, dimension(:,:), allocatable               :: P
        double precision, dimension(:,:), allocatable               :: Q
        double precision, dimension(:,:), allocatable               :: PQ
        double precision, dimension(:,:), allocatable               :: b
        double precision, dimension(:,:), allocatable               :: r
        double precision                                            :: r0
        ! PROCEDURE
        eeen = ones(n,1)
        zeron = zeros(n,1)
        albefa = 0.1d0
        move = 0.5d0

        ! Calculation of the bounds alfa and beta
        zzz1 = low + albefa*(xval-low)
        zzz2 = xval - move*(xmax-xmin)
        zzz  = max(zzz1,zzz2)
        alfa = max(zzz,xmin)
        zzz1 = upp - albefa*(upp-xval)
        zzz2 = xval + move*(xmax-xmin)
        zzz  = min(zzz1,zzz2)
        beta = min(zzz,xmax)

        ! Calculations of p0, q0, P, Q, r and b
        xmami = xmax - xmin
        xmamieps = 0.00001d0*eeen
        xmami = max(xmami,xmamieps)
        xmamiinv = eeen/xmami
        ux1 = upp - xval
        ux2 = ux1*ux1
        xl1 = xval - low
        xl2 = xl1*xl1
        uxinv = eeen/ux1
        xlinv = eeen/xl1
        p0 = zeron
        q0 = zeron
        p0 = max(df0dx,0.0d0)
        q0 = max(-1.0d0*df0dx,0.0d0)
        pq0 = p0 + q0
        p0 = p0 + 0.001d0*pq0
        q0 = q0 + 0.001d0*pq0
        p0 = p0 + raa0*xmamiinv
        q0 = q0 + raa0*xmamiinv
        p0 = p0*ux2
        q0 = q0*xl2
        r0 = f0val - maxval(matmul(transpose(p0),uxinv)) - maxval(matmul(transpose(q0),xlinv))
        !P = zeros(m,n)
        !Q = zeros(m,n)
        P = max(dfdx,0.0d0)
        Q = max(-1.0d0*dfdx,0.0d0)
        PQ = P + Q
        P = P + 0.001d0*PQ
        Q = Q + 0.001d0*PQ
        P = P + matmul(raa,transpose(xmamiinv))
        Q = Q + matmul(raa,transpose(xmamiinv))
        P = MultMatDiag(P,ux2)
        Q = MultMatDiag(Q,xl2)
        r = fval - matmul(P,uxinv) - matmul(Q,xlinv)
        b = -r

        ! Solving the subproblem by a primal-dual Newton method
        call Subsolve(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,m,n,epsimin,low,upp,alfa,beta, &
                    p0,q0,P,Q,a0,a,b,c,d)

        ! Calculations of f0app and fapp.
        ux1 = upp - xmma
        xl1 = xmma - low
        uxinv = eeen/ux1
        xlinv = eeen/xl1
        f0app = r0 + maxval(matmul(transpose(p0),uxinv)) + maxval(matmul(transpose(q0),xlinv))
        fapp  = r + matmul(P,uxinv) + matmul(Q,xlinv)
    end subroutine GCMMAsubDP

    ! Primal-dual Newton method
    subroutine SubSolveReal(xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma, &
                        m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
        ! In/Output Arguments
        real, dimension(:,:), allocatable, intent(out)                :: xmma
        real, dimension(:,:), allocatable, intent(out)                :: ymma
        real, intent(out)                                             :: zmma
        real, dimension(:,:), allocatable, intent(out)                :: lamma
        real, dimension(:,:), allocatable, intent(out)                :: xsimma
        real, dimension(:,:), allocatable, intent(out)                :: etamma
        real, dimension(:,:), allocatable, intent(out)                :: mumma
        real, intent(out)                                             :: zetmma
        real, dimension(:,:), allocatable, intent(out)                :: smma
        integer, intent(in)                                           :: m
        integer, intent(in)                                           :: n
        real, intent(in)                                              :: epsimin
        real, dimension(:,:), allocatable, intent(in)                 :: low
        real, dimension(:,:), allocatable, intent(in)                 :: upp
        real, dimension(:,:), allocatable, intent(in)                 :: alfa
        real, dimension(:,:), allocatable, intent(in)                 :: beta
        real, dimension(:,:), allocatable, intent(in)                 :: p0
        real, dimension(:,:), allocatable, intent(in)                 :: q0
        real, dimension(:,:), allocatable                             :: pq0
        real, dimension(:,:), allocatable, intent(in)                 :: P
        real, dimension(:,:), allocatable, intent(in)                 :: Q
        real, intent(in)                                              :: a0
        real, dimension(:,:), allocatable, intent(in)                 :: a
        real, dimension(:,:), allocatable, intent(in)                 :: b
        real, dimension(:,:), allocatable, intent(in)                 :: c
        real, dimension(:,:), allocatable, intent(in)                 :: d
        ! InternalArguments
        integer                                                       :: itera
        integer                                                       :: ittt
        integer                                                       :: itto
        real                                                          :: epsi
        real                                                          :: residunorm
        real                                                          :: residumax
        real, dimension(:), allocatable                               :: residu1
        real, dimension(:), allocatable                               :: residu2
        real, dimension(:), allocatable                               :: residu
        real, dimension(:,:), allocatable                             :: een
        real, dimension(:,:), allocatable                             :: eem
        real, dimension(:,:), allocatable                             :: epsvecm
        real, dimension(:,:), allocatable                             :: epsvecn
        real, dimension(:,:), allocatable                             :: x
        real, dimension(:,:), allocatable                             :: y
        real                                                          :: z
        real, dimension(:,:), allocatable                             :: dx
        real, dimension(:,:), allocatable                             :: dy
        real                                                          :: dz
        real, dimension(:,:), allocatable                             :: xold
        real, dimension(:,:), allocatable                             :: yold
        real                                                          :: zold
        real, dimension(:,:), allocatable                             :: delx
        real, dimension(:,:), allocatable                             :: dely
        real                                                          :: delz
        real, dimension(:,:), allocatable                             :: lam
        real, dimension(:,:), allocatable                             :: mu
        real, dimension(:,:), allocatable                             :: xsi
        real, dimension(:,:), allocatable                             :: eta
        real, dimension(:,:), allocatable                             :: s
        real, dimension(:,:), allocatable                             :: ux1
        real, dimension(:,:), allocatable                             :: ux2
        real, dimension(:,:), allocatable                             :: ux3
        real, dimension(:,:), allocatable                             :: xl1
        real, dimension(:,:), allocatable                             :: xl2
        real, dimension(:,:), allocatable                             :: xl3
        real, dimension(:,:), allocatable                             :: uxinv1
        real, dimension(:,:), allocatable                             :: xlinv1
        real, dimension(:,:), allocatable                             :: uxinv2
        real, dimension(:,:), allocatable                             :: xlinv2
        real, dimension(:,:), allocatable                             :: dellam
        real, dimension(:,:), allocatable                             :: diagx
        real, dimension(:,:), allocatable                             :: diagy
        real, dimension(:,:), allocatable                             :: diagxinv
        real, dimension(:,:), allocatable                             :: diagyinv
        real, dimension(:,:), allocatable                             :: diaglam
        real, dimension(:,:), allocatable                             :: diaglamyi
        real, dimension(:,:), allocatable                             :: plam
        real, dimension(:,:), allocatable                             :: qlam
        real, dimension(:,:), allocatable                             :: gvec
        real, dimension(:,:), allocatable                             :: dpsidx
        real, dimension(:,:), allocatable                             :: res
        real, dimension(:,:), allocatable                             :: rex
        real, dimension(:,:), allocatable                             :: rey
        real                                                          :: rez
        real, dimension(:,:), allocatable                             :: relam
        real, dimension(:,:), allocatable                             :: rexsi
        real, dimension(:,:), allocatable                             :: reeta
        real, dimension(:,:), allocatable                             :: remu
        real                                                          :: rezet
        real                                                          :: resinew
        real, dimension(:,:), allocatable                             :: blam
        real, dimension(:,:), allocatable                             :: GG
        real, dimension(:,:), allocatable                             :: bb
        real, dimension(:,:), allocatable                             :: bx
        real, dimension(:,:), allocatable                             :: bz
        real, dimension(:,:), allocatable                             :: AA
        real, dimension(:,:), allocatable                             :: Axx
        real, dimension(:,:), allocatable                             :: axz
        real                                                          :: azz
        real, dimension(:,:), allocatable                             :: Alam
        real, dimension(:,:), allocatable                             :: solut
        real, dimension(:,:), allocatable                             :: dlam
        real, dimension(:,:), allocatable                             :: diaglamyiinv
        real, dimension(:,:), allocatable                             :: dellamyi
        real, dimension(:,:), allocatable                             :: dxsi
        real, dimension(:,:), allocatable                             :: deta
        real, dimension(:,:), allocatable                             :: dmu
        real                                                          :: zet
        real                                                          :: dzet
        real, dimension(:,:), allocatable                             :: ds
        real, dimension(:), allocatable                               :: xx
        real, dimension(:), allocatable                               :: dxx
        real, dimension(:), allocatable                               :: stepxx
        real                                                          :: stmxx
        real, dimension(:,:), allocatable                             :: stepalfa
        real                                                          :: stmalfa
        real, dimension(:,:), allocatable                             :: stepbeta
        real                                                          :: stmbeta
        real                                                          :: stmalbe
        real                                                          :: stmalbexx
        real                                                          :: stminv
        real                                                          :: steg
        real, dimension(:,:), allocatable                             :: lamold
        real, dimension(:,:), allocatable                             :: xsiold
        real, dimension(:,:), allocatable                             :: etaold
        real, dimension(:,:), allocatable                             :: muold
        real                                                          :: zetold
        real, dimension(:,:), allocatable                             :: sold
        ! PROCEDURE
        een = ones(n,1)
        eem = ones(m,1)
        epsi = 1.0
        epsvecn = epsi*een
        epsvecm = epsi*eem
        x = 0.5*(alfa+beta)
        y = eem
        z = 1.0
        lam = eem
        xsi = een/(x-alfa)
        xsi = max(xsi,een)
        eta = een/(beta-x)
        eta = max(eta,een)
        mu  = max(eem,0.5*c)
        zet = 1.0
        s = eem
        itera = 0.0
        do while (epsi > epsimin)
            epsvecn = epsi*een
            epsvecm = epsi*eem
            ux1 = upp-x
            xl1 = x-low
            ux2 = ux1*ux1
            xl2 = xl1*xl1
            uxinv1 = een/ux1
            xlinv1 = een/xl1
            plam = p0 + matmul(transpose(P),lam) 
            qlam = q0 + matmul(transpose(Q),lam) 
            gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
            dpsidx = plam/ux2 - qlam/xl2
            rex = dpsidx - xsi + eta
            rey = c + d*y - mu - lam
            rez = a0 - zet - maxval(matmul(transpose(a),lam))
            relam = gvec - a*z - y + s - b
            rexsi = xsi*(x-alfa) - epsvecn
            reeta = eta*(beta-x) - epsvecn
            remu = mu*y - epsvecm
            rezet = zet*z - epsi
            res = lam*s - epsvecm
            residu1 = [rex(:,1),rey(:,1),rez]
            residu2 = [relam(:,1),rexsi(:,1),reeta(:,1),remu(:,1),rezet,res(:,1)]
            residu = [residu1,residu2]
            residunorm = sqrt(dot_product(residu,residu))
            residumax = maxval(abs(residu))
            ittt = 0
            do while ((residumax.gt.0.9*epsi).and.(ittt.lt.200))
                ittt = ittt + 1
                itera = itera + 1
                ux1 = upp - x
                xl1 = x - low
                ux2 = ux1*ux1
                xl2 = xl1*xl1
                ux3 = ux1*ux2
                xl3 = xl1*xl2
                uxinv1 = een/ux1
                xlinv1 = een/xl1
                uxinv2 = een/ux2
                xlinv2 = een/xl2
                plam = p0 + matmul(transpose(P),lam)
                qlam = q0 + matmul(transpose(Q),lam)
                gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
                GG = MultMatDiag(P,uxinv2) - MultMatDiag(Q,xlinv2)
                dpsidx = plam/ux2 - qlam/xl2
                delx = dpsidx - epsvecn/(x-alfa) + epsvecn/(beta-x)
                dely = c + d*y - lam - epsvecm/y
                delz = a0 - maxval(matmul(transpose(a),lam)) - epsi/z
                dellam = gvec - a*z - y - b + epsvecm/lam
                diagx = plam/ux3 + qlam/xl3
                diagx = 2*diagx + xsi/(x-alfa) + eta/(beta-x)
                diagxinv = een/diagx
                diagy = d + mu/y
                diagyinv = eem/diagy
                diaglam = s/lam
                diaglamyi = diaglam + diagyinv
                if (m.lt.n) then
                    blam = dellam + dely/diagy - matmul(GG,(delx/diagx))
                    Alam = diag(diaglamyi) + matmul(MultMatDiag(GG,diagxinv),transpose(GG))
                    ! -- AA Assembly
                    !AA = [Alam     a
                    !       a'   -zet/z] 
                    allocate(AA(size(Alam,1)+1,size(Alam,2)+1))
                    AA((1):(size(Alam,1)),(1):(size(Alam,2))) = Alam
                    AA(1:size(Alam,1),size(Alam,2)+1) = a(:,1)
                    AA(size(Alam,1)+1,1:size(Alam,2)) = a(:,1)
                    AA(size(Alam,1)+1,size(Alam,2)+1) = -zet/z
                    ! -- bb Assembly
                    !bb = [blam' delz]' 
                    allocate(bb(size(blam,1)+1,1))
                    bb(:,1) = [blam(:,1), delz]
                    ! Applying solver
                    solut = Solver(AA,bb)
                    dlam = solut(1:m,:)
                    dz = solut(m+1,1)
                    dx = -delx/diagx - matmul(transpose(GG),dlam)/diagx
                else
                    diaglamyiinv = eem/diaglamyi
                    dellamyi = dellam + dely/diagy
                    Axx = diag(diagx) + matmul(MultMatDiag(GG,diagxinv),transpose(GG))
                    azz = zet/z + maxval(matmul(transpose(a),(a/diaglamyi)))
                    axz = matmul(transpose(-GG),(a/diaglamyi))
                    bx = delx + matmul(transpose(GG),(dellamyi/diaglamyi))
                    bz  = delz - maxval(matmul(transpose(a),(dellamyi/diaglamyi)))
                    ! AA Assembly
                    !AA = [Axx   axz
                    !      axz'  azz] 
                    allocate(AA(size(Axx,1)+1,size(Axx,2)+1))
                    AA((1):(size(Axx,1)),(1):(size(Axx,2))) = Axx
                    AA(1:size(Axx,1),size(Axx,2)+1) = axz(:,1)
                    AA(size(Axx,1)+1,1:size(Axx,2)) = axz(:,1)
                    AA(size(Axx,1)+1,size(Axx,2)+1) = azz
                    ! bb Assembly
                    !bb = [-bx' -bz]' 
                    allocate(bb(size(bx,1)+1,1))
                    bb(:,1) = [-bx(:,1), -bz]
                    ! Applying solver
                    solut = Solver(AA,bb)
                    dx  = solut(1:n,:)
                    dz = solut(n+1,1)
                    dlam = matmul(GG,dx)/diaglamyi - dz*(a/diaglamyi) + dellamyi/diaglamyi
                end if
                deallocate(AA,bb)
                dy = -dely/diagy + dlam/diagy
                dxsi = -xsi + epsvecn/(x-alfa) - (xsi*dx)/(x-alfa)
                deta = -eta + epsvecn/(beta-x) + (eta*dx)/(beta-x)
                dmu  = -mu + epsvecm/y - (mu*dy)/y
                dzet = -zet + epsi/z - zet*dz/z
                ds   = -s + epsvecm/lam - (s*dlam)/lam
                xx  = [y(:,1),z,lam(:,1),xsi(:,1),eta(:,1),mu(:,1),zet,s(:,1)]
                dxx = [dy(:,1),dz,dlam(:,1),dxsi(:,1),deta(:,1),dmu(:,1),dzet,ds(:,1)]
                stepxx = -1.01*dxx/xx
                stmxx  = maxval(stepxx)
                stepalfa = -1.01*dx/(x-alfa)
                stmalfa = maxval(stepalfa)
                stepbeta = 1.01*dx/(beta-x)
                stmbeta = maxval(stepbeta)
                stmalbe  = max(stmalfa,stmbeta)
                stmalbexx = max(stmalbe,stmxx)
                stminv = max(stmalbexx,1.0)
                steg = 1/stminv
                xold = x
                yold = y
                zold = z
                lamold = lam
                xsiold = xsi
                etaold = eta
                muold  = mu
                zetold = zet
                sold = s
                itto = 0
                resinew = 2*residunorm
                do while ((resinew.gt.residunorm).and.(itto.lt.50))
                    itto = itto + 1 
                    x = xold + steg*dx
                    y = yold + steg*dy
                    z = zold + steg*dz
                    lam = lamold + steg*dlam
                    xsi = xsiold + steg*dxsi
                    eta = etaold + steg*deta
                    mu = muold  + steg*dmu
                    zet = zetold + steg*dzet
                    s = sold + steg*ds
                    ux1 = upp - x
                    xl1 = x - low
                    ux2 = ux1*ux1
                    xl2 = xl1*xl1
                    uxinv1 = een/ux1
                    xlinv1 = een/xl1
                    plam = p0 + matmul(transpose(P),lam) 
                    qlam = q0 + matmul(transpose(Q),lam) 
                    gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
                    dpsidx = plam/ux2 - qlam/xl2 
                    rex = dpsidx - xsi + eta
                    rey = c + d*y - mu - lam
                    rez = a0 - zet - maxval(matmul(transpose(a),lam))
                    relam = gvec - a*z - y + s - b
                    rexsi = xsi*(x-alfa) - epsvecn
                    reeta = eta*(beta-x) - epsvecn
                    remu = mu*y - epsvecm
                    rezet = zet*z - epsi
                    res = lam*s - epsvecm
                    residu1 = [rex(:,1), rey(:,1), rez]
                    residu2 = [relam(:,1), rexsi(:,1), reeta(:,1), remu(:,1), rezet, res(:,1)]
                    residu = [residu1, residu2]
                    resinew = sqrt(dot_product(residu,residu))
                    steg = steg/2.0
                end do
                residunorm = resinew
                residumax = maxval(abs(residu))
                steg = 2.0*steg
            end do
            epsi = 0.1*epsi
        end do
        xmma = x
        ymma = y
        zmma = z
        lamma = lam
        xsimma = xsi
        etamma = eta
        mumma = mu
        zetmma = zet
        smma = s
    end subroutine SubSolveReal

    subroutine SubSolveDP(xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma, &
                        m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
        ! In/Output Arguments
        double precision, dimension(:,:), allocatable, intent(out)    :: xmma
        double precision, dimension(:,:), allocatable, intent(out)    :: ymma
        double precision, intent(out)                                 :: zmma
        double precision, dimension(:,:), allocatable, intent(out)    :: lamma
        double precision, dimension(:,:), allocatable, intent(out)    :: xsimma
        double precision, dimension(:,:), allocatable, intent(out)    :: etamma
        double precision, dimension(:,:), allocatable, intent(out)    :: mumma
        double precision, intent(out)                                 :: zetmma
        double precision, dimension(:,:), allocatable, intent(out)    :: smma
        integer, intent(in)                                           :: m
        integer, intent(in)                                           :: n
        double precision, intent(in)                                  :: epsimin
        double precision, dimension(:,:), allocatable, intent(in)     :: low
        double precision, dimension(:,:), allocatable, intent(in)     :: upp
        double precision, dimension(:,:), allocatable, intent(in)     :: alfa
        double precision, dimension(:,:), allocatable, intent(in)     :: beta
        double precision, dimension(:,:), allocatable, intent(in)     :: p0
        double precision, dimension(:,:), allocatable, intent(in)     :: q0
        double precision, dimension(:,:), allocatable                 :: pq0
        double precision, dimension(:,:), allocatable, intent(in)     :: P
        double precision, dimension(:,:), allocatable, intent(in)     :: Q
        double precision, intent(in)                                  :: a0
        double precision, dimension(:,:), allocatable, intent(in)     :: a
        double precision, dimension(:,:), allocatable, intent(in)     :: b
        double precision, dimension(:,:), allocatable, intent(in)     :: c
        double precision, dimension(:,:), allocatable, intent(in)     :: d
        ! InternalArguments
        integer                                                       :: itera
        integer                                                       :: ittt
        integer                                                       :: itto
        double precision                                              :: epsi
        double precision                                              :: residunorm
        double precision                                              :: residumax
        double precision, dimension(:), allocatable                   :: residu1
        double precision, dimension(:), allocatable                   :: residu2
        double precision, dimension(:), allocatable                   :: residu
        double precision, dimension(:,:), allocatable                 :: een
        double precision, dimension(:,:), allocatable                 :: eem
        double precision, dimension(:,:), allocatable                 :: epsvecm
        double precision, dimension(:,:), allocatable                 :: epsvecn
        double precision, dimension(:,:), allocatable                 :: x
        double precision, dimension(:,:), allocatable                 :: y
        double precision                                              :: z
        double precision, dimension(:,:), allocatable                 :: dx
        double precision, dimension(:,:), allocatable                 :: dy
        double precision                                              :: dz
        double precision, dimension(:,:), allocatable                 :: xold
        double precision, dimension(:,:), allocatable                 :: yold
        double precision                                              :: zold
        double precision, dimension(:,:), allocatable                 :: delx
        double precision, dimension(:,:), allocatable                 :: dely
        double precision                                              :: delz
        double precision, dimension(:,:), allocatable                 :: lam
        double precision, dimension(:,:), allocatable                 :: mu
        double precision, dimension(:,:), allocatable                 :: xsi
        double precision, dimension(:,:), allocatable                 :: eta
        double precision, dimension(:,:), allocatable                 :: s
        double precision, dimension(:,:), allocatable                 :: ux1
        double precision, dimension(:,:), allocatable                 :: ux2
        double precision, dimension(:,:), allocatable                 :: ux3
        double precision, dimension(:,:), allocatable                 :: xl1
        double precision, dimension(:,:), allocatable                 :: xl2
        double precision, dimension(:,:), allocatable                 :: xl3
        double precision, dimension(:,:), allocatable                 :: uxinv1
        double precision, dimension(:,:), allocatable                 :: xlinv1
        double precision, dimension(:,:), allocatable                 :: uxinv2
        double precision, dimension(:,:), allocatable                 :: xlinv2
        double precision, dimension(:,:), allocatable                 :: dellam
        double precision, dimension(:,:), allocatable                 :: diagx
        double precision, dimension(:,:), allocatable                 :: diagy
        double precision, dimension(:,:), allocatable                 :: diagxinv
        double precision, dimension(:,:), allocatable                 :: diagyinv
        double precision, dimension(:,:), allocatable                 :: diaglam
        double precision, dimension(:,:), allocatable                 :: diaglamyi
        double precision, dimension(:,:), allocatable                 :: plam
        double precision, dimension(:,:), allocatable                 :: qlam
        double precision, dimension(:,:), allocatable                 :: gvec
        double precision, dimension(:,:), allocatable                 :: dpsidx
        double precision, dimension(:,:), allocatable                 :: res
        double precision, dimension(:,:), allocatable                 :: rex
        double precision, dimension(:,:), allocatable                 :: rey
        double precision                                              :: rez
        double precision, dimension(:,:), allocatable                 :: relam
        double precision, dimension(:,:), allocatable                 :: rexsi
        double precision, dimension(:,:), allocatable                 :: reeta
        double precision, dimension(:,:), allocatable                 :: remu
        double precision                                              :: rezet
        double precision                                              :: resinew
        double precision, dimension(:,:), allocatable                 :: blam
        double precision, dimension(:,:), allocatable                 :: GG
        double precision, dimension(:,:), allocatable                 :: bb
        double precision, dimension(:,:), allocatable                 :: bx
        double precision, dimension(:,:), allocatable                 :: bz
        double precision, dimension(:,:), allocatable                 :: AA
        double precision, dimension(:,:), allocatable                 :: ta
        double precision, dimension(:,:), allocatable                 :: Axx
        double precision, dimension(:,:), allocatable                 :: axz
        double precision                                              :: azz
        double precision, dimension(:,:), allocatable                 :: Alam
        double precision, dimension(:,:), allocatable                 :: solut
        double precision, dimension(:,:), allocatable                 :: dlam
        double precision, dimension(:,:), allocatable                 :: diaglamyiinv
        double precision, dimension(:,:), allocatable                 :: dellamyi
        double precision, dimension(:,:), allocatable                 :: dxsi
        double precision, dimension(:,:), allocatable                 :: deta
        double precision, dimension(:,:), allocatable                 :: dmu
        double precision                                              :: zet
        double precision                                              :: dzet
        double precision, dimension(:,:), allocatable                 :: ds
        double precision, dimension(:), allocatable                   :: xx
        double precision, dimension(:), allocatable                   :: dxx
        double precision, dimension(:), allocatable                   :: stepxx
        double precision                                              :: stmxx
        double precision, dimension(:,:), allocatable                 :: stepalfa
        double precision                                              :: stmalfa
        double precision, dimension(:,:), allocatable                 :: stepbeta
        double precision                                              :: stmbeta
        double precision                                              :: stmalbe
        double precision                                              :: stmalbexx
        double precision                                              :: stminv
        double precision                                              :: steg
        double precision, dimension(:,:), allocatable                 :: lamold
        double precision, dimension(:,:), allocatable                 :: xsiold
        double precision, dimension(:,:), allocatable                 :: etaold
        double precision, dimension(:,:), allocatable                 :: muold
        double precision                                              :: zetold
        double precision, dimension(:,:), allocatable                 :: sold
        ! PROCEDURE
        een = ones(n,1)
        eem = ones(m,1)
        epsi = 1.0d0
        epsvecn = epsi*een
        epsvecm = epsi*eem
        x = 0.5d0*(alfa+beta)
        y = eem
        z = 1
        lam = eem
        xsi = een/(x-alfa)
        xsi = max(xsi,een)
        eta = een/(beta-x)
        eta = max(eta,een)
        mu  = max(eem,0.5*c)
        zet = 1
        s = eem
        itera = 0
        
        do while (epsi > epsimin)
            epsvecn = epsi*een
            epsvecm = epsi*eem
            ux1 = upp-x
            xl1 = x-low
            ux2 = ux1*ux1
            xl2 = xl1*xl1
            uxinv1 = een/ux1
            xlinv1 = een/xl1
            plam = p0 + matmul(transpose(P),lam) 
            qlam = q0 + matmul(transpose(Q),lam) 
            gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
            dpsidx = plam/ux2 - qlam/xl2
            rex = dpsidx - xsi + eta
            rey = c + d*y - mu - lam
            rez = a0 - zet - maxval(matmul(transpose(a),lam))
            relam = gvec - a*z - y + s - b
            rexsi = xsi*(x-alfa) - epsvecn
            reeta = eta*(beta-x) - epsvecn
            remu = mu*y - epsvecm
            rezet = zet*z - epsi
            res = lam*s - epsvecm
            residu1 = [rex(:,1),rey(:,1),rez]
            residu2 = [relam(:,1),rexsi(:,1),reeta(:,1),remu(:,1),rezet,res(:,1)]
            residu = [residu1,residu2]
            residunorm = sqrt(dot_product(residu,residu))
            residumax = maxval(abs(residu))
            ittt = 0

            do while ((residumax.gt.0.9d0*epsi).and.(ittt.lt.200))
                ittt = ittt + 1
                itera = itera + 1
                ux1 = upp - x
                xl1 = x - low
                ux2 = ux1*ux1
                xl2 = xl1*xl1
                ux3 = ux1*ux2
                xl3 = xl1*xl2
                uxinv1 = een/ux1
                xlinv1 = een/xl1
                uxinv2 = een/ux2
                xlinv2 = een/xl2
                plam = p0 + matmul(transpose(P),lam)
                qlam = q0 + matmul(transpose(Q),lam)
                gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
                GG = MultMatDiag(P,uxinv2) - MultMatDiag(Q,xlinv2)
                dpsidx = plam/ux2 - qlam/xl2
                delx = dpsidx - epsvecn/(x-alfa) + epsvecn/(beta-x)
                dely = c + d*y - lam - epsvecm/y
                delz = a0 - maxval(matmul(transpose(a),lam)) - epsi/z
                dellam = gvec - a*z - y - b + epsvecm/lam
                diagx = plam/ux3 + qlam/xl3
                diagx = 2*diagx + xsi/(x-alfa) + eta/(beta-x)
                diagxinv = een/diagx
                diagy = d + mu/y
                diagyinv = eem/diagy
                diaglam = s/lam
                diaglamyi = diaglam + diagyinv

                if (m.lt.n) then
                    blam = dellam + dely/diagy - matmul(GG,(delx/diagx))
                    Alam = diag(diaglamyi) + matmul(MultMatDiag(GG,diagxinv),transpose(GG))
                    ! -- AA Assembly
                    !AA = [Alam     a
                    !       a'   -zet/z] 
                    allocate(AA(size(Alam,1)+1,size(Alam,2)+1))
                    AA((1):(size(Alam,1)),(1):(size(Alam,2))) = Alam
                    AA(1:size(Alam,1),size(Alam,2)+1) = a(:,1)
                    AA(size(Alam,1)+1,1:size(Alam,2)) = a(:,1)
                    AA(size(Alam,1)+1,size(Alam,2)+1) = -zet/z
                    ! -- bb Assembly
                    !bb = [blam' delz]' 
                    allocate(bb(size(blam,1)+1,1))
                    bb(:,1) = [blam(:,1), delz]
                    ! Applying solver
                    solut = Solver(AA,bb)
                    dlam = solut(1:m,:)
                    dz = solut(m+1,1)
                    dx = -delx/diagx - matmul(transpose(GG),dlam)/diagx
                else
                    diaglamyiinv = eem/diaglamyi
                    dellamyi = dellam + dely/diagy
                    Axx = diag(diagx) + matmul(MultMatDiag(GG,diagxinv),transpose(GG))
                    azz = zet/z + maxval(matmul(transpose(a),(a/diaglamyi)))
                    axz = matmul(transpose(-GG),(a/diaglamyi))
                    bx = delx + matmul(transpose(GG),(dellamyi/diaglamyi))
                    bz  = delz - maxval(matmul(transpose(a),(dellamyi/diaglamyi)))
                    ! AA Assembly
                    !AA = [Axx   axz
                    !      axz'  azz] 
                    allocate(AA(size(Axx,1)+1,size(Axx,2)+1))
                    AA((1):(size(Axx,1)),(1):(size(Axx,2))) = Axx
                    AA(1:size(Axx,1),size(Axx,2)+1) = axz(:,1)
                    AA(size(Axx,1)+1,1:size(Axx,2)) = axz(:,1)
                    AA(size(Axx,1)+1,size(Axx,2)+1) = azz
                    ! bb Assembly
                    !bb = [-bx' -bz]' 
                    allocate(bb(size(bx,1)+1,1))
                    bb(:,1) = [-bx(:,1), -bz]
                    ! Applying solver
                    solut = Solver(AA,bb)
                    dx  = solut(1:n,:)
                    dz = solut(n+1,1)
                    dlam = matmul(GG,dx)/diaglamyi - dz*(a/diaglamyi) + dellamyi/diaglamyi
                end if
                deallocate(AA,bb)

                dy = -dely/diagy + dlam/diagy
                dxsi = -xsi + epsvecn/(x-alfa) - (xsi*dx)/(x-alfa)
                deta = -eta + epsvecn/(beta-x) + (eta*dx)/(beta-x)
                dmu  = -mu + epsvecm/y - (mu*dy)/y
                dzet = -zet + epsi/z - zet*dz/z
                ds   = -s + epsvecm/lam - (s*dlam)/lam

                xx  = [y(:,1),z,lam(:,1),xsi(:,1),eta(:,1),mu(:,1),zet,s(:,1)]
                dxx = [dy(:,1),dz,dlam(:,1),dxsi(:,1),deta(:,1),dmu(:,1),dzet,ds(:,1)]

                stepxx = -1.01d0*dxx/xx
                stmxx  = maxval(stepxx)
                stepalfa = -1.01d0*dx/(x-alfa)
                stmalfa = maxval(stepalfa)
                stepbeta = 1.01d0*dx/(beta-x)
                stmbeta = maxval(stepbeta)
                stmalbe  = max(stmalfa,stmbeta)
                stmalbexx = max(stmalbe,stmxx)
                stminv = max(stmalbexx,1.0d0)
                steg = 1/stminv

                xold = x
                yold = y
                zold = z
                lamold = lam
                xsiold = xsi
                etaold = eta
                muold  = mu
                zetold = zet
                sold = s

                itto = 0
                resinew = 2*residunorm

                do while ((resinew.gt.residunorm).and.(itto.lt.50))
                    itto = itto + 1 
                    x = xold + steg*dx
                    y = yold + steg*dy
                    z = zold + steg*dz
                    lam = lamold + steg*dlam
                    xsi = xsiold + steg*dxsi
                    eta = etaold + steg*deta
                    mu = muold  + steg*dmu
                    zet = zetold + steg*dzet
                    s = sold + steg*ds
                    ux1 = upp - x
                    xl1 = x - low
                    ux2 = ux1*ux1
                    xl2 = xl1*xl1
                    uxinv1 = een/ux1
                    xlinv1 = een/xl1
                    plam = p0 + matmul(transpose(P),lam) 
                    qlam = q0 + matmul(transpose(Q),lam) 
                    gvec = matmul(P,uxinv1) + matmul(Q,xlinv1)
                    dpsidx = plam/ux2 - qlam/xl2 
                    rex = dpsidx - xsi + eta
                    rey = c + d*y - mu - lam
                    rez = a0 - zet - maxval(matmul(transpose(a),lam))
                    relam = gvec - a*z - y + s - b
                    rexsi = xsi*(x-alfa) - epsvecn
                    reeta = eta*(beta-x) - epsvecn
                    remu = mu*y - epsvecm
                    rezet = zet*z - epsi
                    res = lam*s - epsvecm
                    residu1 = [rex(:,1), rey(:,1), rez]
                    residu2 = [relam(:,1), rexsi(:,1), reeta(:,1), remu(:,1), rezet, res(:,1)]
                    residu = [residu1, residu2]
                    resinew = sqrt(dot_product(residu,residu))
                    steg = steg/2.0d0
                end do
                residunorm = resinew
                residumax = maxval(abs(residu))
                steg = 2.0d0*steg
            end do
            epsi = 0.1d0*epsi
        end do
        xmma = x
        ymma = y
        zmma = z
        lamma = lam
        xsimma = xsi
        etamma = eta
        mumma = mu
        zetmma = zet
        smma = s
    end subroutine SubSolveDP

    ! KKT Checkk
    subroutine KKTCheckReal(residu,residunorm,residumax,x,y,z,lam,xsi,eta,mu,zet,s, &
                        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        implicit none
        ! In/Output Arguments
        real, dimension(:,:), allocatable, intent(inout)  :: x
        real, dimension(:,:), allocatable, intent(inout)  :: y
        real, intent(inout)                               :: z
        real, dimension(:,:), allocatable, intent(inout)  :: lam
        real, dimension(:,:), allocatable, intent(inout)  :: xsi
        real, dimension(:,:), allocatable, intent(inout)  :: eta
        real, dimension(:,:), allocatable, intent(inout)  :: mu
        real, intent(inout)                               :: zet
        real, dimension(:,:), allocatable, intent(inout)  :: s
        real, dimension(:,:), allocatable, intent(inout)  :: xmin
        real, dimension(:,:), allocatable, intent(inout)  :: xmax
        real, dimension(:,:), allocatable, intent(inout)  :: df0dx
        real, dimension(:,:), allocatable, intent(inout)  :: fval
        real, dimension(:,:), allocatable, intent(inout)  :: dfdx
        real, intent(inout)                               :: a0
        real, dimension(:,:), allocatable, intent(inout)  :: a
        real, dimension(:,:), allocatable, intent(inout)  :: c
        real, dimension(:,:), allocatable, intent(inout)  :: d
        real, dimension(:), allocatable, intent(out)      :: residu
        real, intent(out)                                 :: residunorm
        real, intent(out)                                 :: residumax
        ! internal variables
        real, dimension(:,:), allocatable                 :: rex
        real, dimension(:,:), allocatable                 :: rey
        real                                              :: rez
        real, dimension(:), allocatable                   :: residu1
        real, dimension(:), allocatable                   :: residu2
        real, dimension(:,:), allocatable                 :: relam
        real, dimension(:,:), allocatable                 :: rexsi
        real, dimension(:,:), allocatable                 :: reeta
        real, dimension(:,:), allocatable                 :: remu
        real                                              :: rezet
        real, dimension(:,:), allocatable                 :: res
        ! process
        rex   = df0dx + maxval(matmul(transpose(dfdx),lam)) - xsi + eta
        rey   = c + d*y - mu - lam
        rez   = a0 - zet - maxval(matmul(transpose(a),lam))
        relam = fval - a*z - y + s

        rexsi = xsi*(x - xmin)
        reeta = eta*(xmax - x)
        remu  = mu*y
        rezet = zet*z
        res   = lam*s
        ! results
        residu1 = [rex(:,1), rey(:,1), rez]
        residu2 = [relam(:,1), rexsi(:,1), reeta(:,1), remu(:,1), rezet, res(:,1)]
        residu = [residu1, residu2]
        residunorm = sqrt(dot_product(residu,residu)) 
        residumax = maxval(abs(residu)) 
    end subroutine KKTCheckReal

    subroutine KKTCheckDP(residu,residunorm,residumax,x,y,z,lam,xsi,eta,mu,zet,s, &
                    xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        implicit none
        ! In/Output Arguments                                                                   
        double precision, dimension(:,:), allocatable, intent(inout)  :: x
        double precision, dimension(:,:), allocatable, intent(inout)  :: y
        double precision, intent(inout)                               :: z
        double precision, dimension(:,:), allocatable, intent(inout)  :: lam
        double precision, dimension(:,:), allocatable, intent(inout)  :: xsi
        double precision, dimension(:,:), allocatable, intent(inout)  :: eta
        double precision, dimension(:,:), allocatable, intent(inout)  :: mu
        double precision, intent(inout)                               :: zet
        double precision, dimension(:,:), allocatable, intent(inout)  :: s
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmin
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmax
        double precision, dimension(:,:), allocatable, intent(inout)  :: df0dx
        double precision, dimension(:,:), allocatable, intent(inout)  :: fval
        double precision, dimension(:,:), allocatable, intent(inout)  :: dfdx
        double precision, intent(inout)                               :: a0
        double precision, dimension(:,:), allocatable, intent(inout)  :: a
        double precision, dimension(:,:), allocatable, intent(inout)  :: c
        double precision, dimension(:,:), allocatable, intent(inout)  :: d
        double precision, dimension(:), allocatable, intent(out)      :: residu
        double precision, intent(out)                                 :: residunorm
        double precision, intent(out)                                 :: residumax
        ! internal variables
        double precision, dimension(:,:), allocatable                 :: rex
        double precision, dimension(:,:), allocatable                 :: rey
        double precision                                              :: rez
        double precision, dimension(:), allocatable                   :: residu1
        double precision, dimension(:), allocatable                   :: residu2
        double precision, dimension(:,:), allocatable                 :: relam
        double precision, dimension(:,:), allocatable                 :: rexsi
        double precision, dimension(:,:), allocatable                 :: reeta
        double precision, dimension(:,:), allocatable                 :: remu
        double precision                                              :: rezet
        double precision, dimension(:,:), allocatable                 :: res
        ! process
        rex   = df0dx + maxval(matmul(transpose(dfdx),lam)) - xsi + eta
        rey   = c + d*y - mu - lam
        rez   = a0 - zet - maxval(matmul(transpose(a),lam))
        relam = fval - a*z - y + s
        rexsi = xsi*(x - xmin)
        reeta = eta*(xmax - x)
        remu  = mu*y
        rezet = zet*z
        res   = lam*s
        ! results
        residu1 = [rex(:,1), rey(:,1), rez]
        residu2 = [relam(:,1), rexsi(:,1), reeta(:,1), remu(:,1), rezet, res(:,1)]
        residu = [residu1, residu2]
        residunorm = sqrt(dot_product(residu,residu)) 
        residumax = maxval(abs(residu)) 
    end subroutine KKTCheckDP

    ! Concheck
    subroutine ConcheckReal(conserv,m,epsimin,f0app,f0valnew,fapp,fvalnew)
        implicit none
        ! In/Output Arguments
        integer, intent(inout)                              :: m
        integer, intent(inout)                              :: conserv
        real, intent(inout)                                 :: epsimin
        real, intent(inout)                                 :: f0valnew
        real, intent(inout)                                 :: f0app
        real, dimension(:,:), allocatable, intent(inout)    :: fapp
        real, dimension(:,:), allocatable, intent(inout)    :: fvalnew
        ! internal variables
        real                                                :: f0appe
        real, dimension(:,:), allocatable                   :: fappe
        real, dimension(:,:), allocatable                   :: eeem
        conserv = 0
        eeem = ones(m,1)
        f0appe = f0app + epsimin
        fappe = fapp + epsimin*eeem
        if (all(fappe.ge.fvalnew).and.(f0appe.ge.f0valnew)) then
            conserv = 1
        end if
    end subroutine ConcheckReal

    subroutine ConcheckDP(conserv,m,epsimin,f0app,f0valnew,fapp,fvalnew)
        implicit none
        ! In/Output Arguments
        integer, intent(inout)                                        :: m
        integer, intent(inout)                                        :: conserv
        double precision, intent(inout)                               :: epsimin
        double precision, intent(inout)                               :: f0valnew
        double precision, intent(inout)                               :: f0app
        double precision, dimension(:,:), allocatable, intent(inout)  :: fapp
        double precision, dimension(:,:), allocatable, intent(inout)  :: fvalnew
        ! internal variables
        double precision                                              :: f0appe
        double precision, dimension(:,:), allocatable                 :: fappe
        double precision, dimension(:,:), allocatable                 :: eeem
        conserv = 0
        eeem = ones(m,1)
        f0appe = f0app + epsimin
        fappe = fapp + epsimin*eeem
        if (all(fappe.ge.fvalnew).and.(f0appe.ge.f0valnew)) then
            conserv = 1
        end if
    end subroutine ConcheckDP

    ! Asymp
    subroutine AsympReal(low,upp,raa0,raa,outeriter,n,xval,xold1,xold2,xmin,xmax, &
                        raa0eps,raaeps,df0dx,dfdx)
        implicit none
        ! In/Output Arguments
        integer, intent(inout)                            :: n
        integer, intent(inout)                            :: outeriter
        real, intent(inout)                               :: raa0
        real, intent(inout)                               :: raa0eps
        real, dimension(:,:), allocatable, intent(inout)  :: raa
        real, dimension(:,:), allocatable, intent(inout)  :: raaeps
        real, dimension(:,:), allocatable, intent(inout)  :: low
        real, dimension(:,:), allocatable, intent(inout)  :: upp
        real, dimension(:,:), allocatable, intent(inout)  :: xval
        real, dimension(:,:), allocatable, intent(inout)  :: xold1
        real, dimension(:,:), allocatable, intent(inout)  :: xold2
        real, dimension(:,:), allocatable, intent(inout)  :: xmin
        real, dimension(:,:), allocatable, intent(inout)  :: xmax
        real, dimension(:,:), allocatable, intent(inout)  :: dfdx
        real, dimension(:,:), allocatable, intent(inout)  :: df0dx
        ! internal variables
        real                                              :: asyinit
        real                                              :: asyincr
        real                                              :: asydecr
        real, dimension(:,:), allocatable                 :: xmami
        real, dimension(:,:), allocatable                 :: lowmin
        real, dimension(:,:), allocatable                 :: lowmax
        real, dimension(:,:), allocatable                 :: uppmin
        real, dimension(:,:), allocatable                 :: uppmax
        real, dimension(:,:), allocatable                 :: zzz
        real, dimension(:,:), allocatable                 :: factor
        real, dimension(:,:), allocatable                 :: eeen
        real, dimension(:,:), allocatable                 :: xmamieps
        eeen = ones(n,1)
        xmami = xmax - xmin 
        xmamieps = 0.00001*eeen
        xmami = max(xmami,xmamieps)
        raa0 = maxval(matmul(transpose(abs(df0dx)),xmami))
        raa0 = maxval([raa0eps,(0.1/n)*raa0])
        raa = matmul(abs(dfdx),xmami)
        raa = max(raaeps,(0.1/n)*raa)
        asyinit = 0.5       
        asyincr = 1.2       
        asydecr = 0.7
        ! Calculation of the asymptotes low and upp
        if (outeriter.lt.2.5) then
            low = xval - asyinit*(xmami)
            upp = xval + asyinit*(xmami)
        else
            zzz = (xval-xold1)*(xold1-xold2)
            factor = eeen
            where (zzz.gt.0.0)
                factor = asyincr
            elsewhere (zzz.lt.0.0)
                factor = asydecr
            end where
            low = xval - factor*(xold1 - low)
            upp = xval + factor*(upp - xold1)
            lowmin = xval - 10.0*(xmami)
            lowmax = xval - 0.01*(xmami)
            uppmin = xval + 0.01*(xmami)
            uppmax = xval + 10.0*(xmami)
            low = max(low,lowmin)
            low = min(low,lowmax)
            upp = min(upp,uppmax)
            upp = max(upp,uppmin)
        end if
    end subroutine AsympReal

    subroutine AsympDP(low,upp,raa0,raa,outeriter,n,xval,xold1,xold2,xmin,xmax, &
                    raa0eps,raaeps,df0dx,dfdx)
        implicit none
        ! In/Output Arguments
        integer, intent(inout)                                        :: n
        integer, intent(inout)                                        :: outeriter
        double precision, intent(inout)                               :: raa0
        double precision, intent(inout)                               :: raa0eps
        double precision, dimension(:,:), allocatable, intent(inout)  :: raa
        double precision, dimension(:,:), allocatable, intent(inout)  :: raaeps
        double precision, dimension(:,:), allocatable, intent(inout)  :: low
        double precision, dimension(:,:), allocatable, intent(inout)  :: upp
        double precision, dimension(:,:), allocatable, intent(inout)  :: xval
        double precision, dimension(:,:), allocatable, intent(inout)  :: xold1
        double precision, dimension(:,:), allocatable, intent(inout)  :: xold2
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmin
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmax
        double precision, dimension(:,:), allocatable, intent(inout)  :: dfdx
        double precision, dimension(:,:), allocatable, intent(inout)  :: df0dx
        ! internal variables
        double precision                                              :: asyinit
        double precision                                              :: asyincr
        double precision                                              :: asydecr
        double precision, dimension(:,:), allocatable                 :: xmami
        double precision, dimension(:,:), allocatable                 :: lowmin
        double precision, dimension(:,:), allocatable                 :: lowmax
        double precision, dimension(:,:), allocatable                 :: uppmin
        double precision, dimension(:,:), allocatable                 :: uppmax
        double precision, dimension(:,:), allocatable                 :: zzz
        double precision, dimension(:,:), allocatable                 :: factor
        double precision, dimension(:,:), allocatable                 :: eeen
        double precision, dimension(:,:), allocatable                 :: xmamieps
        eeen = ones(n,1)
        xmami = xmax - xmin 
        xmamieps = 0.00001d0*eeen
        xmami = max(xmami,xmamieps)
        raa0 = maxval(matmul(transpose(abs(df0dx)),xmami))
        raa0 = maxval([raa0eps,(0.1d0/n)*raa0])
        raa = matmul(abs(dfdx),xmami)
        raa = max(raaeps,(0.1d0/n)*raa)
        asyinit = 0.5d0
        asyincr = 1.2d0
        asydecr = 0.7d0
        ! Calculation of the asymptotes low and upp
        if (outeriter.lt.2.5) then
            low = xval - asyinit*(xmami)
            upp = xval + asyinit*(xmami)
        else
            zzz = (xval-xold1)*(xold1-xold2)
            factor = eeen
            where (zzz.gt.0.0d0)
                factor = asyincr
            elsewhere (zzz.lt.0.0d0)
                factor = asydecr
            end where
            low = xval - factor*(xold1 - low)
            upp = xval + factor*(upp - xold1)
            lowmin = xval - 10.0d0*(xmami)
            lowmax = xval - 0.01d0*(xmami)
            uppmin = xval + 0.01d0*(xmami)
            uppmax = xval + 10.0d0*(xmami)
            low = max(low,lowmin)
            low = min(low,lowmax)
            upp = min(upp,uppmax)
            upp = max(upp,uppmin)
        end if
    end subroutine AsympDP

    ! raaUpdate
    subroutine raaUpdateReal(raa0,raa,xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, &
                            f0app,fapp,raa0eps,raaeps,epsimin)
        implicit none
        ! In/Output Arguments
        real, intent(inout)                                          :: raa0
        real, intent(inout)                                          :: raa0eps
        real, intent(inout)                                          :: f0app
        real, intent(inout)                                          :: f0valnew
        real, intent(inout)                                          :: epsimin
        real, dimension(:,:), allocatable, intent(inout)             :: raa
        real, dimension(:,:), allocatable, intent(inout)             :: raaeps
        real, dimension(:,:), allocatable, intent(inout)             :: xmma
        real, dimension(:,:), allocatable, intent(inout)             :: xval
        real, dimension(:,:), allocatable, intent(inout)             :: fapp
        real, dimension(:,:), allocatable, intent(inout)             :: fvalnew
        real, dimension(:,:), allocatable, intent(inout)             :: low
        real, dimension(:,:), allocatable, intent(inout)             :: upp
        real, dimension(:,:), allocatable, intent(inout)             :: xmin
        real, dimension(:,:), allocatable, intent(inout)             :: xmax
        ! internal variables
        integer                                                      :: m
        integer                                                      :: n
        real                                                         :: raacofmin
        real                                                         :: f0appe
        real                                                         :: zz0
        real                                                         :: deltaraa0
        real                                                         :: raacof
        real, dimension(:,:), allocatable                            :: zzz
        real, dimension(:,:), allocatable                            :: eeem
        real, dimension(:,:), allocatable                            :: eeen
        real, dimension(:,:), allocatable                            :: xmami
        real, dimension(:,:), allocatable                            :: xmamieps
        real, dimension(:,:), allocatable                            :: xxux
        real, dimension(:,:), allocatable                            :: xxxl
        real, dimension(:,:), allocatable                            :: xxul
        real, dimension(:,:), allocatable                            :: ulxx
        real, dimension(:,:), allocatable                            :: fappe
        real, dimension(:,:), allocatable                            :: fdelta
        real, dimension(:,:), allocatable                            :: deltaraa
        ! process
        raacofmin = 1.0E-12
        m = size(raa)
        n = size(xmma)
        eeem = ones(m,1)
        eeen = ones(n,1)
        xmami = xmax - xmin 
        xmamieps = 0.00001*eeen
        xmami = max(xmami,xmamieps)
        !
        xxux = (xmma - xval)/(upp - xmma)
        xxxl = (xmma - xval)/(xmma - low)
        xxul = xxux*xxxl
        ulxx = (upp - low)/xmami
        raacof = maxval([matmul(transpose(xxul),ulxx),raacofmin])
        !
        f0appe = f0app + 0.5*epsimin 
        if (f0valnew.gt.f0appe) then
            deltaraa0 = (1.0/raacof)*(f0valnew - f0app)
            zz0 = 1.1*(raa0 + deltaraa0)
            zz0 = min(zz0,10.0*raa0)
            !zz0 = min(zz0,1000.0*raa0)
            raa0 = zz0
        end if 
        !
        fappe = fapp + 0.5*epsimin*eeem
        fdelta = fvalnew - fappe
        deltaraa = (1/raacof)*(fvalnew - fapp)
        zzz = 1.1*(raa + deltaraa)
        zzz = min(zzz,10.0*raa)
        !zzz = min(zzz,1000.0*raa)
        where (fdelta.gt.0.0)
            raa = zzz
        end where
    end subroutine raaUpdateReal

    subroutine raaUpdateDP(raa0,raa,xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, &
                            f0app,fapp,raa0eps,raaeps,epsimin)
        implicit none
        ! In/Output Arguments
        double precision, intent(inout)                               :: raa0
        double precision, intent(inout)                               :: raa0eps
        double precision, intent(inout)                               :: f0app
        double precision, intent(inout)                               :: f0valnew
        double precision, intent(inout)                               :: epsimin
        double precision, dimension(:,:), allocatable, intent(inout)  :: raa
        double precision, dimension(:,:), allocatable, intent(inout)  :: raaeps
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmma
        double precision, dimension(:,:), allocatable, intent(inout)  :: xval
        double precision, dimension(:,:), allocatable, intent(inout)  :: fapp
        double precision, dimension(:,:), allocatable, intent(inout)  :: fvalnew
        double precision, dimension(:,:), allocatable, intent(inout)  :: low
        double precision, dimension(:,:), allocatable, intent(inout)  :: upp
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmin
        double precision, dimension(:,:), allocatable, intent(inout)  :: xmax
        ! internal variables
        integer                                                       :: m
        integer                                                       :: n
        double precision                                              :: raacofmin
        double precision                                              :: f0appe
        double precision                                              :: zz0
        double precision                                              :: deltaraa0
        double precision                                              :: raacof
        double precision, dimension(:,:), allocatable                 :: zzz
        double precision, dimension(:,:), allocatable                 :: eeem
        double precision, dimension(:,:), allocatable                 :: eeen
        double precision, dimension(:,:), allocatable                 :: xmami
        double precision, dimension(:,:), allocatable                 :: xmamieps
        double precision, dimension(:,:), allocatable                 :: xxux
        double precision, dimension(:,:), allocatable                 :: xxxl
        double precision, dimension(:,:), allocatable                 :: xxul
        double precision, dimension(:,:), allocatable                 :: ulxx
        double precision, dimension(:,:), allocatable                 :: fappe
        double precision, dimension(:,:), allocatable                 :: fdelta
        double precision, dimension(:,:), allocatable                 :: deltaraa
        ! process
        raacofmin = 1.0d-12
        m = size(raa)
        n = size(xmma)
        eeem = ones(m,1)
        eeen = ones(n,1)
        xmami = xmax - xmin
        xmamieps = 0.00001d0*eeen
        xmami = max(xmami,xmamieps)
        !
        xxux = (xmma - xval)/(upp - xmma)
        xxxl = (xmma - xval)/(xmma - low)
        xxul = xxux*xxxl
        ulxx = (upp - low)/xmami
        raacof = maxval([matmul(transpose(xxul),ulxx),raacofmin])
        !
        f0appe = f0app + 0.5d0*epsimin
        if (f0valnew.gt.f0appe) then
            deltaraa0 = (1.0d0/raacof)*(f0valnew - f0app)
            zz0 = 1.1d0*(raa0 + deltaraa0)
            zz0 = min(zz0,10.0d0*raa0)
            !zz0 = min(zz0,1000.0d0*raa0)
            raa0 = zz0
        end if 
        !
        fappe = fapp + 0.5d0*epsimin*eeem
        fdelta = fvalnew - fappe
        deltaraa = (1.0d0/raacof)*(fvalnew - fapp)
        zzz = 1.1d0*(raa + deltaraa)
        zzz = min(zzz,10.0d0*raa)
        !zzz = min(zzz,1000.0d0*raa)
        where (fdelta.gt.0.0d0)
            raa = zzz
        end where
    end subroutine raaUpdateDP

    ! Solver for [A]{x}={B}
    function SolverReal(A,B) result(X)
        implicit none
        ! In/Output Arguments
        real, dimension(:,:), allocatable, intent(in)       :: B
        real, dimension(:,:), allocatable, intent(in)       :: A
        ! internal variables
        integer                                             :: n,m
        integer                                             :: info
        integer, dimension(:), allocatable                  :: ipiv
        real, dimension(:,:), allocatable                   :: X
        n = size(B,1)
        m = size(B,2)
        allocate(ipiv(n))
        call SGESV(n,m,A,n,ipiv,B,n,info)
        if(info.ne.0) then
            write(unit=*, fmt=*) "ERROR in DGESV solver, failed with info", info 
            stop "Check the system"
        end if
        X = B
    end function SolverReal

    function SolverDP(A,B) result(X)
        implicit none
        ! In/Output Arguments
        double precision, dimension(:,:), allocatable, intent(in)   :: B
        double precision, dimension(:,:), allocatable, intent(in)   :: A
        ! internal variables
        integer                                                     :: n,m
        integer                                                     :: info
        integer, dimension(:), allocatable                          :: ipiv
        double precision, dimension(:,:), allocatable               :: X
        n = size(B,1)
        m = size(B,2)
        allocate(ipiv(n))
        call DGESV(n,m,A,n,ipiv,B,n,info)
        if(info.ne.0) then
            write(unit=*, fmt=*) "ERROR in DGESV solver, failed with info", info 
            stop "Check the system"
        end if
        X = B
    end function SolverDP
end module GCMMA_Routines

