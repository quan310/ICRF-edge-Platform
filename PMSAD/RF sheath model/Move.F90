Module Move2D
    Use TypeModule
    Implicit none
    contains	
	subroutine moveB(InputParticle,NX,NY,Ex,Ey,Bx,By,Bz,AccFactor,AccFactorB)
		Implicit none
		Integer i,N1,N2
		Type(ParticleBundle),intent(inout) :: InputParticle		
		Integer,intent(in) :: NX,NY
		real(8),intent(in) :: Ex(NX,NY),Ey(NX,NY),Bx(NX,NY),By(NX,NY),Bz(NX,NY)
		real(8) :: ExTemp(NX,NY),EyTemp(NX,NY),BxTemp(NX,NY),ByTemp(NX,NY),BzTemp(NX,NY)
		real(8) :: AccFactor,AccFactorB
		real(8) :: X1,X2,Y1,Y2,E1,E2,B1,B2,B3,F11,F12,F21,F22,Vx1,Vy1,Vz1,s
		REAL(8) :: TT(3)		
		ExTemp = Ex*AccFactor*2.d0
		EyTemp = Ey*AccFactor*2.d0
		BxTemp = Bx*AccFactorB
		ByTemp = By*AccFactorB
		BzTemp = Bz*AccFactorB		
		do i = 1,InputParticle%Npar
			N1 = ceiling(InputParticle%PhaseSpace(i)%X)
			X1 = DBLE(N1) - InputParticle%PhaseSpace(i)%X
			X2 = 1.d0-X1			
			N2 = ceiling(InputParticle%PhaseSpace(i)%Y)
			Y1 = DBLE(N2) - InputParticle%PhaseSpace(i)%Y
			Y2 = 1.d0-Y1			
			F11 = X1*Y1
			F12 = X1*Y2
			F21 = X2*Y1
			F22 = X2*Y2		
			E1 = ExTemp(N1,N2)*F11 + ExTemp(N1+1,N2)*F21 + ExTemp(N1,N2+1)*F12 + ExTemp(N1+1,N2+1)*F22			
			E2 = EyTemp(N1,N2)*F11 + EyTemp(N1+1,N2)*F21 + EyTemp(N1,N2+1)*F12 + EyTemp(N1+1,N2+1)*F22
			InputParticle%PhaseSpace(i)%VX = InputParticle%PhaseSpace(i)%Vx + E1/2.0
			InputParticle%PhaseSpace(i)%VY = InputParticle%PhaseSpace(i)%VY + E2/2.0
			B1 = BxTemp(N1,N2)*F11 + BxTemp(N1+1,N2)*F21 + BxTemp(N1,N2+1)*F12 + BxTemp(N1+1,N2+1)*F22
			B2 = ByTemp(N1,N2)*F11 + ByTemp(N1+1,N2)*F21 + ByTemp(N1,N2+1)*F12 + ByTemp(N1+1,N2+1)*F22
			B3 = BzTemp(N1,N2)*F11 + BzTemp(N1+1,N2)*F21 + BzTemp(N1,N2+1)*F12 + BzTemp(N1+1,N2+1)*F22
			tt(1) = B1
			tt(2) = B2
			tt(3) = B3
			Vx1 = InputParticle%PhaseSpace(i)%VX - InputParticle%PhaseSpace(i)%VZ*tt(2) + InputParticle%PhaseSpace(i)%Vy*tt(3)
			Vy1 = InputParticle%PhaseSpace(i)%VY + InputParticle%PhaseSpace(i)%VZ*tt(1) - InputParticle%PhaseSpace(i)%Vx*tt(3)
			Vz1 = InputParticle%PhaseSpace(i)%VZ + InputParticle%PhaseSpace(i)%VX*tt(2) - InputParticle%PhaseSpace(i)%VY*tt(1)		
			s = tt(1)**2+tt(2)**2+tt(3)**2		
			InputParticle%PhaseSpace(i)%VX = InputParticle%PhaseSpace(i)%VX-Vz1*2.0*tt(2)/(1.0+s)+Vy1*2.0*tt(3)/(1.0+s)
			InputParticle%PhaseSpace(i)%VY = InputParticle%PhaseSpace(i)%VY+Vz1*2.0*tt(1)/(1.0+s)-Vx1*2.0*tt(3)/(1.0+s)
			InputParticle%PhaseSpace(i)%VZ = InputParticle%PhaseSpace(i)%VZ+Vx1*2.0*tt(2)/(1.0+s)-Vy1*2.0*tt(1)/(1.0+s)
			InputParticle%PhaseSpace(i)%VX = InputParticle%PhaseSpace(i)%VX+E1/2.0
			InputParticle%PhaseSpace(i)%VY = InputParticle%PhaseSpace(i)%VY+E2/2.0
			InputParticle%PhaseSpace(i)%X = InputParticle%PhaseSpace(i)%X+InputParticle%PhaseSpace(i)%VX
			InputParticle%PhaseSpace(i)%Y = InputParticle%PhaseSpace(i)%Y+InputParticle%PhaseSpace(i)%VY	
		end do
		return
	end subroutine moveB			                              
End Module Move2D