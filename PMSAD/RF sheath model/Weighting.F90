 Module WeightingModule
    Use TypeModule 
	Use ParticleModule
    implicit none
    contains 
	subroutine Weighting2D(InputParticle,Nx,Ny,Rho)
                implicit none
                Type(ParticleBundle),intent(in) :: InputParticle
                Integer(4),intent(in) :: Nx,Ny
                real(8),intent(out) :: Rho(Nx,Ny)
                real(8) :: X1,X2,Y1,Y2
                Integer(4) :: i,N1,N2
                Rho=0.d0
                do i=1,InputParticle%Npar
                   N1=Ceiling(InputParticle%PhaseSpace(i)%X)
                   X1=Dble(N1)-InputParticle%PhaseSpace(i)%X
                   X2=1.d0-X1
                   N2=Ceiling(InputParticle%PhaseSpace(i)%Y) 
                   Y1=Dble(N2)-InputParticle%PhaseSpace(i)%Y
                   Y2=1.d0-Y1
if(N1>Nx.or.N1<1.or.N2>Ny.or.N2<1) then
write(*,*) 'N1,N2',N1,N2
else
                   Rho(N1,N2)=Rho(N1,N2)+X1*Y1
                   Rho(N1,N2+1)=Rho(N1,N2+1)+X1*Y2
                   Rho(N1+1,N2)=Rho(N1+1,N2)+X2*Y1
                   Rho(N1+1,N2+1)=Rho(N1+1,N2+1)+X2*Y2
end if
               end do
  Rho(1:Nx,1)=2.d0*Rho(1:Nx,1) 
  Rho(1:Nx,Ny)=2.d0*Rho(1:Nx,Ny) 
             return
  end subroutine Weighting2D
        
    Subroutine Sencondaries2(InputParticle,Sez,Ser)
            Implicit none
            Type(ParticleBundle),intent(inout) :: InputParticle      
            Type(ParticleOne) :: ParticleTemp
            Integer(4) :: i
            Real(8) :: Sez,Ser  
            ParticleTemp%X=Sez
            ParticleTemp%Y=Ser
            Call MaxwellianSe3(ParticleTemp)
            ParticleTemp%Vz=ParticleTemp%Vz*1.d0/InputParticle%VFactor
            ParticleTemp%Vx=ParticleTemp%Vx*1.d0/InputParticle%VFactor
            ParticleTemp%Vy=ParticleTemp%Vy*1.d0/InputParticle%VFactor
            Call AddParticle(ParticleTemp,InputParticle)
            return 
    End Subroutine Sencondaries2

subroutine MaxwellianSe3(ParticleTemp)
               implicit none              
               Type(ParticleOne),intent(inout) :: ParticleTemp
               real(8) ::  Mass=9.1095d-31,Temprature=3.d0*11605.d0
               real(8) :: V,Beta,FuncA,FuncB
               real(8) :: Theta,CosTheta,SinTheta,Fai               
               Beta=1.d0/(1.3807d-23*Temprature)
               FuncA=1.d0
               FuncB=0.d0
              do while(FuncA>FuncB)
                 Call DRandom(R)
                      FuncA=R*R
                 Call DRandom(R)
                      FuncB=-exp*R*Dlog(R)
              end do
              V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)		 
              ParticleTemp%Vx=0.5*V
              ParticleTemp%Vy=0.5*V
              ParticleTemp%Vz=0.5*V
             return 
           end subroutine MaxwellianSe3  

subroutine againDione(InputParticleion,InputParticlee,Nx,Ny,rhoagain)
                implicit none
                Type(ParticleBundle),intent(inout) :: InputParticleion,InputParticlee
                Integer(4),intent(in) :: Nx,Ny
                real(8) :: X1,X2,Y1,Y2,againx,againy
                Integer(4) :: i,N1,N2,rhoagain(1:Ny),k,Nchou
                Real(8) :: Rnum,Rnum1,Rnum2,Rnum3     
            do k=1,Ny			   
                do i=1,rhoagain(k)                  
				    Call DRandom(Rnum)
                    if(Rnum<0.9) then                        
                         Call DRandom(Rnum1)
                         Call DRandom(Rnum3)
                          Nchou=int(Rnum1*InputParticlee%Npar)
                          if(Nchou<=InputParticlee%Npar.and.Nchou>=1) then
                               if(InputParticlee%PhaseSpace(Nchou)%X>Nx-3)   then
                                 14    if(Nchou<InputParticlee%Npar)   then
                                                   Nchou=Nchou+1
                                                   if(InputParticlee%PhaseSpace(Nchou)%X>Nx-3) go to 14
                                          end if
                                          if(InputParticlee%PhaseSpace(Nchou)%X>Nx-3) then
                                                Nchou=int(Rnum1*InputParticlee%Npar)
                                               13    if(Nchou>1)   then
                                                             Nchou=Nchou-1
                                                            if(InputParticlee%PhaseSpace(Nchou)%X>Nx-3) go to 13
                                                       end if
                                         end if
                                  end if
                               end if
                          if(InputParticlee%PhaseSpace(Nchou)%X>Nx-3) then
                               write(*,*) 'shaixuan'
                          end if
                               if(InputParticlee%PhaseSpace(Nchou)%X>=Nx-2) then                      
                                      againx=InputParticlee%PhaseSpace(Nchou)%X-Rnum3
                                else  
                                      againx=InputParticlee%PhaseSpace(Nchou)%X+Rnum3
                                end if
                               if(InputParticlee%PhaseSpace(Nchou)%X>=Nx-1) then
                                   againx=Nx-1.001
                               end if
                               Call DRandom(Rnum2)
                               againy=k-1+Rnum2
                               call Sencondaries2(InputParticlee,againx,againy) 
                               call Sencondaries2(InputParticleion,againx,againy) 
                   end if
                 end do
              end do
             return
  end subroutine againDione           
                           
End Module WeightingModule