Module InitilalizationModule
	Use FileIO
	Use Parallel
	Use Constants
	Implicit none
    contains
	
	Subroutine  ParticleInit(InputGas,NSpecy,InputSpecy,OutputPICParticlePara,OutputParticleBundle,dx,dt,Density,Nx,Ny,PPG,NProcess)  
            Implicit none
            Type(Gas), intent(in) :: InputGas 
            Integer(4), intent(in) :: NSpecy
            Type(OneSpecy), intent(in) ::  InputSpecy(1:NSpecy)
            Type(PICParticlePara),intent(out) :: OutputPICParticlePara(1:NSpecy)
            Type(ParticleBundle),intent(out) :: OutputParticleBundle(1:NSpecy)
            Integer(4),intent(in) :: Nx,Ny,PPG,NProcess
            Real(8),intent(in) :: dx,dt,Density
            Logical :: Status
            Integer(4) :: i
				
				do i=1,NSpecy
                    Call PICParticleParaInit(InputSpecy(i),OutputPICParticlePara(i),dx,dt,Density,PPG)
                    OutputParticleBundle(i)%XFactor=dx
                    OutputParticleBundle(i)%VFactor=dx/dt
                    OutputParticleBundle(i)%Mass=InputSpecy(i)%Mass
                    OutputParticleBundle(i)%Charge=InputSpecy(i)%Charge  
                    If (NProcess==1) then
                         Call LoadParticle(InputSpecy(i),OutputParticleBundle(i),Status)
                    Else                   
                         Call ParallelLoadParticle(InputSpecy(i),OutputParticleBundle(i),Status)
                    End If
                    If (Status) then
                         Call PhaseSpaceInit2D(InputGas,InputSpecy(i),OutputParticleBundle(i),Nx,Ny,PPG,Nprocess,Density,dx)
                    else
                         Call UpdatePositionBundle(OutputParticleBundle(i),-1_4)
                    End If
                    Call UpdateVelocityBundle(OutputParticleBundle(i),-1_4) 
                End do
            return
        End Subroutine  ParticleInit
       
             Subroutine PICParticleParaInit(InputSpecy,OutputPICParticlePara,dx,dt,Density,PPG)
                Implicit none
                Type(OneSpecy), intent(in) :: InputSpecy
                Type(PICParticlePara),intent(out) :: OutputPICParticlePara
                Real(8), intent(in) :: dx,dt,Density
                Integer(4) :: PPG,SubCycle
                Real(8) :: Charge,QdM
                SubCycle=1
                OutputPICParticlePara%SubCycle=SubCycle
                OutputPICParticlePara%Timer=1
                 
                Charge=dble(InputSpecy%Charge)*ElectronCharge
                QdM=Charge/InputSpecy%Mass
                OutputPICParticlePara%RhoFactor=Charge*Density/dble(PPG)  
                OutputPICParticlePara%Wq=Density/dble(PPG)
                OutputPICParticlePara%AccFactor=0.5d0*QdM*dt/(dx/dt)
                OutputPICParticlePara%AccFactorB=0.5d0*QdM*dt
                OutputPICParticlePara%QFactor=Charge
                OutputPICParticlePara%Value=0.d0
                OutputPICParticlePara%Temp=0.d0 
               return
             End Subroutine  PICParticleParaInit
          
           
             Subroutine PhaseSpaceInit2D(InputGas,InputSpecy,OutputParticleBundle,Nx,Ny,PPG,Nprocss,Density,dx)
                      Use MaxwellModule
                      Implicit none
                      Type(Gas), intent(in) :: InputGas
                      Type(OneSpecy),intent(in) :: InputSpecy 
                      Type(ParticleBundle),intent(out) :: OutputParticleBundle
                      Integer(4),intent(in) :: Nx,Ny,PPG,Nprocss
                      Real(8),intent(in) ::  Density,dx
                      Type(ParticleOne) :: ParticleTemp
                      Integer(4) :: i,j,k,NParticle,PPGLocal
                      Real(8) :: XMin,XMax,YMin,YMax
                      OutputParticleBundle%Name=InputSpecy%Name
                      OutputParticleBundle%NPar=0
                      PPGLocal=PPG/Nprocss 
                      do i=1,Nx-1
                           do j=1,Ny-1
                                 XMin=Dble(i-1)
                                 XMax=XMin+1.d0
                                 YMin=Dble(j-1)  !-0.5d0
                                 YMax=YMin+1.d0
                                 If (YMin<0.d0) YMin=0.d0
                                  do k=1,PPGLocal
                                        Call RandomPosition(ParticleTemp%X,XMin,XMax)
                                        Call RandomPosition(ParticleTemp%Y,YMin,YMax)
                                        Call Maxwellian(InputGas,ParticleTemp) 
                                        Call AddParticle(ParticleTemp,OutputParticleBundle)
                                 end do
                          end do
                      end do      
                      return
                      contains
                          Subroutine RandomPosition(X,MinPosition,MaxPosition)
                              Implicit none
                              Real(8),intent(in) ::  MaxPosition,MinPosition 
                              Real(8),intent(out) ::  X
                              Call DRandom(R)
                              X=MinPosition+R*(MaxPosition-MinPosition)
                              return 
                           end  Subroutine RandomPosition 
                    end subroutine PhaseSpaceInit2D

subroutine add(InputGas,leftup,ParticleGlobal,Xmin,Xmax,Ymin,Ymax,dx,nn)		
		Use MaxwellModule
		Use Input
        Implicit none
		
		Type(Gas), intent(in) :: InputGas
		Integer(4) :: k,i,j
		
		Integer(4),intent(in) ::  leftup(1:(InputNx),1:(InputNy)),nn
		Type(ParticleBundle),intent(inout) :: ParticleGlobal
		Real(8),intent(in):: XMin,XMax,YMin,YMax,dx									
		Type(ParticleOne) :: ParticleTemp	
		Real(8) :: VFactor,XFactor
		
		do i = 1,inputNx-1
	              do j= 1,inputNy-1
				if (leftup(i,j)>0) then
					
					do k = 1,leftup(i,j)
		
						Call DRandom(R)
ParticleTemp%X=R+i-1
Call DRandom(R)
ParticleTemp%Y=R+j-1
if(ParticleTemp%X>Xmax) then
write(*,*) '错误',ParticleTemp%X,ParticleTemp%Y
end if
						Call Maxwellian2(InputGas,ParticleTemp) 
						
						VFactor = 1.0/ParticleGlobal%VFactor
						ParticleTemp%Vx = VFactor*ParticleTemp%Vx
						ParticleTemp%Vy = VFactor*ParticleTemp%Vy
						ParticleTemp%Vz = VFactor*ParticleTemp%Vz

						Call AddParticle(ParticleTemp,ParticleGlobal)
					end do
				end if
		end do
		end do

		return
      		
	end subroutine add

End Module InitilalizationModule

