       subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
		
		double precision xJ1,xJ2,xJ3,xJ, xB(3,3), U(3,3), xB2(3,3),UT(3,3)
		double precision sigma(3,3),xIdentity(3,3),lamda1_2,lamda2_2,lamda3_2
		
		xMu = props(1)
		
		sigma=0.
		
		xIdentity=0.0
		xIdentity(1,1)=1.0
		xIdentity(2,2)=1.0
		xIdentity(3,3)=1.0
		
		xB=0.0
		xB2=0.0

		! loop through all blocks
		do km = 1, nblock

		U(1,1) = stretchNew(km,1)
		U(2,2) = stretchNew(km,2)
		U(3,3) = stretchNew(km,3)
		U(1,2) = stretchNew(km,4)
		U(2,1) = U(1,2)
		
		if (nshr .eq. 1) then
			U(2,3) = 0.0
			U(1,3) = 0.0
		else
			U(2,3) = stretchNew(km,5)
			U(1,3) = stretchNew(km,6)
		end if
		
		U(3,2)=U(2,3)
		U(3,1)=U(1,3)
		
		! calculate J
		xJ1 = U(1,1)*(U(2,2)*U(3,3)-U(2,3)**2)
		xJ2 = -U(1,2)*(U(2,1)*U(3,3)-U(2,3)*U(3,1)) 
		xJ3 = U(1,3)*(U(2,1)*U(3,2)-U(2,2)*U(3,1)) 
		xJ= xJ1+xJ2+xJ3	

		! Define the square of the deviatoric stretch tensor (B)
		CALL KMLT(U,U,xB)
		!xB=matmul(U,U)
		
		! Define the forth power of the deviatoric stretch tensor (B * B)
		CALL KMLT(xB,xB,xB2)
		!xB2=matmul(xB,xB)
				
		! Define some constants to simplify the stress expression
		lamda1_2=xB(1,1)
		lamda2_2=xB(2,2)
		lamda3_2=xB(3,3)
		
		! I-1
		I1 = (lamda1_2+lamda2_2+lamda3_2)
		
		! I-2
		I2 =( lamda1_2*lamda2_2 + lamda2_2*lamda3_2 + lamda3_2*lamda1_2 )  
		
		!new stress matrix
		sigma = (xMu/xJ**3) * (I1*xB - xB2 - (I2-xJ**3)*xIdentity)
		
		!The corotational stress
		stressNew(km,1) = sigma(1,1)
		stressNew(km,2) = sigma(2,2)
		stressNew(km,3) = sigma(3,3)
		stressNew(km,4) = sigma(1,2)
		if (nshr .eq. 3) then
			stressNew(km,5) = sigma(2,3)
			stressNew(km,6) = sigma(3,1)
		end if
		
		end do
		
      return
      end
	  

	  
	  
	  SUBROUTINE KMLT(DM1,DM2,DM)
C      
      include 'vaba_param.inc'
	  
      PARAMETER (M=3,N=3)
      double precision DM1(M,N),DM2(M,N),DM(M,N),X
C
      DO 10 I=1,M
      DO 10 J=1,N
      X=0.0
      DO 20 K=1,M 
      X=X+DM1(I,K)*DM2(K,J)
20    CONTINUE
      DM(I,J)=X
10    CONTINUE
      RETURN
      END