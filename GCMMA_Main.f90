program Example
    use GCMMA_variables
    use GCMMA_Routines
    implicit none
    ! --------------------------------------------------------------- !
    ! these subroutines are essentially an adaptation/translation of  !
    ! the code presented by Professor Svanberg originally programmed  !
    ! in Matlab. The definition of the model is presented below and   !
    ! in the subroutine GCMMA_variables.f90 (avoids modifying other   !
    ! elements of the code).                                          !
    !                                                                 !
    ! Note: for the execution of the code it is necessary to use the  !
    !       Lapack library to execute the SGESV and DGESV routines    !
    !       and additionally a makefile is presented to facilitate    !
    !       the compilation and execution.                            !
    ! --------------------------------------------------------------- !
    m = 2
    n = 3
    allocate(xval(n,1),xmin(n,1),xmax(n,1))
    allocate(a(m,1),c(m,1),d(m,1),raa(m,1),raaeps(m,1))
    allocate(df0dx(n,1),fval(m,1),fvalnew(m,1),dfdx(m,n))
    epsimin = 0.0000001
    outeriter = 0
    maxoutit = 1
    kkttol = 0
    xval(:,1) = [4.0,3.0,2.0]
    xold1 = xval
    xold2 = xval
    xmin(:,1) = [0.0,0.0,0.0]
    xmax(:,1) = [5.0,5.0,5.0]
    low = xmin
    upp = xmax
    a0 = 1.0
    raa0 = 0.01
    raa0eps = 0.000001
    a(:,1) = [0.0,0.0]
    d(:,1) = [1.0,1.0]
    c(:,1) = [1000.0,1000.0]
    raa(:,1) = 0.01*[1, 1]
    raaeps(:,1)  = 0.000001*[1, 1]
    ! --------------------------------------------------------------- !
    !               Here star the GCMMA optimization process          !
    ! --------------------------------------------------------------- !
    if (outeriter.lt.0.5) then
        call ObjectiveFunction(xval,f0val,df0dx,fval,dfdx)
        innerit = 0
        outvector1 = [real(outeriter),real(innerit), xval]
        outvector2 = [f0val, fval]
    end if
    ! The outer iterations start
    kktnorm = kkttol + 10
    outit = 0        
    do while ((kktnorm.gt.kkttol).and.(outit.lt.maxoutit))
        ! counter
        outit = outit + 1
        outeriter = outeriter + 1
        ! The parameters low, upp, raa0 and raa are calculated
        call Asymp(low,upp,raa0,raa,outeriter,n,xval,xold1,xold2,xmin,xmax, &
                    raa0eps,raaeps,df0dx,dfdx)
        ! The GCMMA subproblem is solved at the point xval
        call GCMMA_Sub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp,m,n, &
                    outeriter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val, &
                    df0dx,fval,dfdx,a0,a,c,d)
        call ObjectiveFunction(xmma,f0valnew,df0dx,fvalnew,dfdx)
        ! It is checked if the approximations are conservative
        call Concheck(conserv,m,epsimin,f0app,f0valnew,fapp,fvalnew)
        
        ! While the approximations are non-conservative (conserv=0)
        innerit = 0
        if (conserv.eq.0) then
            do while ((conserv.eq.0).and.(innerit.le.15))
                innerit = innerit + 1
                ! New values on the parameters raa0 and raa are calculated
                call raaUpdate(raa0,raa,xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, &
                                f0app,fapp,raa0eps,raaeps,epsimin)
                ! The GCMMA subproblem is solved with these new raa0 and raa
                call GCMMA_Sub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp,m,n, &
                            outeriter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val, &
                            df0dx,fval,dfdx,a0,a,c,d)
                call ObjectiveFunction(xmma,f0valnew,df0dx,fvalnew,dfdx)
                ! It is checked if the approximations have become conservative
                call Concheck(conserv,m,epsimin,f0app,f0valnew,fapp,fvalnew)
            end do
        end if
        ! Updating of some vectors
        xold2 = xold1
        xold1 = xval
        xval  = xmma
        ! Getting the f0val, df0dx, fval and dfdx
        call ObjectiveFunction(xval,f0val,df0dx,fval,dfdx)
        ! Residual vector of KKT conditions
        call kktcheck(residu,kktnorm,residumax,xmma,ymma,zmma,lam,xsi,eta, &
                        mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        outvector1 = [real(outeriter),real(innerit), xval(:,1)]
        outvector2 = [f0val, fval(:,1)]
    end do
    ! Print results
    write(unit=*, fmt=*) 'outvector1'
    write(unit=*, fmt=*) outvector1
    write(unit=*, fmt=*) 'outvector2'
    write(unit=*, fmt=*) outvector2
    write(unit=*, fmt=*) 'xval'
    write(unit=*, fmt=*) xval
end program Example
