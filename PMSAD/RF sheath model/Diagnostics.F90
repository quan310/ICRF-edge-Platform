Module Diagnostics
    Use TypeModule
    Use ParticleModule 
    Use Parallel 
    Use FileIO 
    Implicit none
    contains

    Subroutine  DensityTemperatureInit(NSpecy,InputParticleBundle,InputNx,InputNy,PICParaGlobal)
         Use Parallel,only : Nprocess
         Implicit none
         Integer(4),intent(in) ::  NSpecy,InputNx,InputNy
         Type(ParticleBundle),intent(in) ::  InputParticleBundle(NSpecy)
         Type(PICParticlePara),intent(in) ::  PICParaGlobal(NSpecy)
        
         Integer(4),parameter :: NSpecyMax=10_4 
         Type(Grid),save :: DensityGrid(NSpecyMax),TemperatureGrid(NSpecyMax),TempGrid,TempGrid2
         Integer(4),save :: i,j,Nx,Ny,Nxy,Timer
         Real(8),save  :: Mass(NSpecyMax),VFactor(NSpecyMax)
            
            Timer=0  
            Nx=InputNx
            Ny=InputNy
            Nxy=Nx*Ny  
            TempGrid%Nx=Nx
            TempGrid%Ny=Ny
            TempGrid2%Nx=Nx
            TempGrid2%Ny=Ny
         do i=1,NSpecy
              Write(DensityGrid(i)%Name,*) trim(InputParticleBundle(i)%Name),'De'
              Write(TemperatureGrid(i)%Name,*) trim(InputParticleBundle(i)%Name),'Te'
              VFactor(i)=InputParticleBundle(i)%VFactor 
              DensityGrid(i)%Nx=Nx
              DensityGrid(i)%Ny=Ny
              DensityGrid(i)%Value=0.d0
              TemperatureGrid(i)%Nx=Nx
              TemperatureGrid(i)%Ny=Ny
              TemperatureGrid(i)%Value=0.d0
         End do     
         Return          
    Entry  DensityTemperature(NSpecy,InputParticleBundle)
        Timer=Timer+1  
        Do i=1, NSpecy
              Call WeightingEnergy(InputParticleBundle(i),Nx,Ny,TempGrid%Value,TempGrid2%Value)
              DensityGrid(i)%Value(1:Nxy)=DensityGrid(i)%Value(1:Nxy)+TempGrid%Value(1:Nxy)
              TemperatureGrid(i)%Value(1:Nxy)=TemperatureGrid(i)%Value(1:Nxy)+TempGrid2%Value(1:Nxy)  
        End do
        Return    
    Entry  DensityTemperatureFinal(NSpecy,PICParaGlobal)
          Do i=1,NSpecy 
                Call MPI_Reduce(DensityGrid(i)%Value,TempGrid%Value,Nxy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)
                Call MPI_Reduce(TemperatureGrid(i)%Value,TempGrid2%Value,Nxy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)
                !Write(*,*)  'ssa',TempGrid%Value(3000),TempGrid2%Value(3000),Timer
                TempGrid%Name=DensityGrid(i)%Name
                TempGrid2%Name=TemperatureGrid(i)%Name
                TempGrid%Value=TempGrid%Value/Dble(Timer)
                TempGrid2%Value=VFactor(i)*VFactor(i)*TempGrid2%Value/Dble(Timer)/Nprocess
                do j=1,Nxy
                     IF(TempGrid%Value(j)>0.d0) then
                     TempGrid2%Value(j)=TempGrid2%Value(j)/TempGrid%Value(j)
                     End if
                End do
                TempGrid%Value=TempGrid%Value*PICParaGlobal(i)%wQ
                
                If(Myid==Master) Call DumpGrid(TempGrid)
                If(Myid==Master) Call DumpGrid(TempGrid2)
          End do
        Return 
    End Subroutine  DensityTemperatureInit 

     subroutine WeightingEnergy(InputParticle,Nx,Ny,Rho,Temperature)
                implicit none
                Type(ParticleBundle),intent(in) :: InputParticle
                
                Integer(4),intent(in) :: Nx,Ny
                real(8),intent(out) :: Rho(Nx,Ny),Temperature(Nx,Ny)
                real(8) :: X1,X2,Y1,Y2,Energy
                Integer(4) :: i,j,N1,N2
                Rho=0.d0
                Temperature=0.d0
                
                !Write(*,*) InputParticle%Mass
                do i=1,InputParticle%Npar
                   Call CalEnergy(InputParticle%Mass,InputParticle%PhaseSpace(i),Energy)
                   !Write(*,*) Energy
                   N1=Ceiling(InputParticle%PhaseSpace(i)%X)
                   X1=Dble(N1)-InputParticle%PhaseSpace(i)%X
                   X2=1.d0-X1
                   N2=Ceiling(InputParticle%PhaseSpace(i)%Y) 
                   Y1=Dble(N2)-InputParticle%PhaseSpace(i)%Y
                   Y2=1.d0-Y1
                   Rho(N1,N2)=Rho(N1,N2)+X1*Y1
                   Rho(N1,N2+1)=Rho(N1,N2+1)+X1*Y2
                   Rho(N1+1,N2)=Rho(N1+1,N2)+X2*Y1
                   Rho(N1+1,N2+1)=Rho(N1+1,N2+1)+X2*Y2
                  
                   Temperature(N1,N2)=Temperature(N1,N2)+X1*Y1*Energy
                   Temperature(N1,N2+1)=Temperature(N1,N2+1)+X1*Y2*Energy
                   Temperature(N1+1,N2)=Temperature(N1+1,N2)+X2*Y1*Energy
                   Temperature(N1+1,N2+1)=Temperature(N1+1,N2+1)+X2*Y2*Energy
               end do
               Rho(1:Nx,1)=2.d0*Rho(1:Nx,1) 
                Temperature(1:Nx,1)=2.d0*Temperature(1:Nx,1) 
               
               
             return
  end subroutine WeightingEnergy
  
       subroutine AvaPhiEInit(Phi,Ex,Ey)
                implicit none
                Type(Grid),intent(in) :: Phi,Ex,Ey
                Type(Grid),save :: TempGrid,TempGrid2,TempGrid3,TempGridEy,TempGridEx
                Integer(4),save :: Nx,Ny,Nxy,Timer
                Nx=Phi%Nx
                Ny=Phi%Ny
                Nxy=Nx*Ny
                Timer=0
                Write(TempGrid%Name,*) trim(Phi%Name),'Ava'                
                TempGrid%Nx=Phi%Nx
                TempGrid%Ny=Phi%Ny
                TempGrid%Value(1:Nxy)=0.d0
                Write(TempGrid2%Name,*) trim(Ex%Name),'Ava'  
                TempGrid2%Nx=Ex%Nx
                TempGrid2%Ny=Ex%Ny
                TempGrid2%Value(1:Nxy)=0.d0
                Write(TempGrid3%Name,*) trim(Ey%Name),'Ava'                
                TempGrid3%Nx=Ey%Nx
                TempGrid3%Ny=Ey%Ny
                TempGrid3%Value(1:Nxy)=0.d0
                Write(TempGridEx%Name,*) 'EXst'          
                TempGridEX%Nx=Ex%Nx
                TempGridEX%Ny=Ex%Ny
                TempGridEx%Value(1:Nxy)=0.d0

             Write(TempGridEy%Name,*) 'Eyst'          
                TempGridEy%Nx=Ey%Nx
                TempGridEy%Ny=Ey%Ny
                TempGridEy%Value(1:Nxy)=0.d0
                
                return
        Entry AvaPhiE(Phi,Ex,Ey)
                  Timer=Timer+1
                  TempGrid%Value(1:Nxy)=TempGrid%Value(1:Nxy)+Phi%Value(1:Nxy)
                  TempGrid2%Value(1:Nxy)=TempGrid2%Value(1:Nxy)+Ex%Value(1:Nxy)
                  TempGrid3%Value(1:Nxy)=TempGrid3%Value(1:Nxy)+Ey%Value(1:Nxy)
                  TempGridEy%Value(1:Nxy)=Ey%Value(1:Nxy)
                  TempGridEx%Value(1:Nxy)=Ex%Value(1:Nxy)

       return  
        Entry AvaPhiEFinal()
                TempGrid%Value=TempGrid%Value/Dble(Timer) 
                TempGrid2%Value=TempGrid2%Value/Dble(Timer) 
                TempGrid3%Value=TempGrid3%Value/Dble(Timer)
                If(Myid==Master) Call DumpGrid(TempGrid)
                If(Myid==Master) Call DumpGrid(TempGrid2)
                If(Myid==Master) Call DumpGrid(TempGrid3)
                return   
       End subroutine AvaPhiEInit





   Subroutine  CalPowerDensityInit(NSpecy,InputParticleBundle,InputNx,InputNy,Inputdx,PICParaGlobal)
         USE input,only:Inputdt
         USE Constants,only:Epsilon
         Use ElectronModule
         Use HGasModule
         implicit none
         Integer(4),intent(in) ::  InputNx,InputNy,NSpecy
         Type(ParticleBundle),intent(in) ::  InputParticleBundle(0:NSpecy)
         Type(PICParticlePara),intent(in) ::  PICParaGlobal(0:NSpecy)
         Type(Grid),intent(in) :: Ex,Ey
         Real(8),intent(in) :: Inputdx
         Integer(4),parameter :: NSpecyMax=2_4 
         Type(Grid),save :: PowerDensity,TempGrid,TempGrid2,TempGrid3
         Integer(4),save :: i,j,Nx,Ny,Nxy,Timer
         Real(8),save  :: dx,VFactor(0:NSpecyMax)
         Type(Grid),save :: powergGrid(0:NSpecyMax)
         Real(8):: PL,PH
         Type(Grid),save :: Ex1Grid,ExGrid,Ey1Grid,powerelec,powerion,ExtotalGrid,powert,powerte,powertion,ExGridall
            dx=Inputdx
            Timer=0
            Nx=InputNx
            Ny=InputNy
            Nxy=Nx*Ny
            PowerDensity%Name='PowerDensity' 
            PowerDensity%Nx=Nx
            PowerDensity%Ny=Ny
            PowerDensity%Value=0.d0
            
            TempGrid%Nx=Nx
            TempGrid%Ny=Ny
            TempGrid%Value=0.d0
            
            TempGrid2%Name='PowerDensity'    
            TempGrid2%Nx=Nx
            TempGrid2%Ny=Ny
            TempGrid2%Value=0.d0

           TempGrid3%Name='PowerDensity'    
            TempGrid3%Nx=Nx
            TempGrid3%Ny=Ny
            TempGrid3%Value=0.d0

            Ex1Grid%Nx=Nx
            Ex1Grid%Ny=Ny
            Ex1Grid%Value(1:Nxy)=0.d0

            ExGrid%Name='Exstdisp'   
            ExGrid%Nx=Nx
            ExGrid%Ny=Ny
            ExGrid%Value(1:Nxy)=0.d0

            ExGridall%Name='Exstdisp'   
            ExGridall%Nx=Nx
            ExGridall%Ny=Ny
            ExGridall%Value(1:Nxy)=0.d0

            ExtotalGrid%Name='displacement' 
            ExtotalGrid%Nx=Nx
            ExtotalGrid%Ny=Ny
            ExtotalGrid%Value(1:Nxy)=0.d0
            
            Ey1Grid%Nx=Nx
            Ey1Grid%Ny=Ny
            Ey1Grid%Value(1:Nxy)=0.d0
          
           powerelec%Name='powerelec' 
            powerelec%Nx=Nx
            powerelec%Ny=Ny
            powerelec%Value(1:Nxy)=0.d0

            powerion%Name='powerion' 
            powerion%Nx=Nx
            powerion%Ny=Ny
            powerion%Value(1:Nxy)=0.d0

            powert%Name='powert' 
            powert%Nx=Nx
            powert%Ny=Ny
            powert%Value(1:Nxy)=0.d0

            powerte%Name='powerte' 
            powerte%Nx=Nx
            powerte%Ny=Ny
            powerte%Value(1:Nxy)=0.d0

            powertion%Name='powertion' 
            powertion%Nx=Nx
            powertion%Ny=Ny
            powertion%Value(1:Nxy)=0.d0

            do i=0,NSpecy
              Write(powergGrid(i)%Name,*) 'DiPower',trim(InputParticleBundle(i)%Name)

              VFactor(i)=InputParticleBundle(i)%VFactor 
              
              powergGrid(i)%Nx=Nx
              powergGrid(i)%Ny=Ny
              powergGrid(i)%Value=0.d0
             end do

           

            Return  
      Entry  CalPowerDensity(NSpecy,InputParticleBundle,Ex,Ey,PICParaGlobal)
        Timer=Timer+1  

        do i=0,NSpecy
         Call WeightingPowerDensity(InputParticleBundle(i),Nx,Ny,Ex%Value,Ey%Value,TempGrid%Value,dx)
  
               If(i==0) then
                       PowerDensity%Value=PowerDensity%Value+TempGrid%Value*VFactor(i)*PICParaGlobal(i)%wQ
                       powergGrid(i)%Value=powergGrid(i)%Value+TempGrid%Value*VFactor(i)*PICParaGlobal(i)%wQ
                         powerelec%Value=powerelec%Value+TempGrid%Value*VFactor(i)*PICParaGlobal(i)%wQ
                        TempGrid3%Value=TempGrid%Value*VFactor(i)*PICParaGlobal(i)%wQ

               end if
                            
                           PowerDensity%Value=PowerDensity%Value+TempGrid%Value*dble(ArSpecy(i)%Charge)*VFactor(i)*PICParaGlobal(i)%wQ
                           powergGrid(i)%Value=powergGrid(i)%Value+TempGrid%Value*dble(ArSpecy(i)%Charge)*VFactor(i)*PICParaGlobal(i)%wQ
                           powerion%Value=powerion%Value+TempGrid%Value*dble(ArSpecy(i)%Charge)*VFactor(i)*PICParaGlobal(i)%wQ

                          TempGrid3%Value=TempGrid%Value*dble(ArSpecy(i)%Charge)*VFactor(i)*PICParaGlobal(i)%wQ

        end do


      do i=1,Nxy
     
        ExGrid%Value(i)=(Ex%Value(i)-Ex1Grid%Value(i))/Inputdt*Epsilon*Ex%Value(i)+(Ey%Value(i)-Ey1Grid%Value(i))/Inputdt*Epsilon*Ey%Value(i)
     end do
  
       ! call Elec(ExGrid%Value,TempGrid%Value,Ny,Nx,EpsilonGlobal%Value)
         ExtotalGrid%Value=ExtotalGrid%Value+ExGrid%Value
        Ex1Grid%Value=Ex%Value
        Ey1Grid%Value=Ey%Value

        Return
      Entry  CalPowerDensityFinal( )

           Call MPI_Reduce(powerelec%Value,TempGrid%Value,Nxy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)   
    
         TempGrid%Value=TempGrid%Value/Dble(Timer) 
         TempGrid%name=powerelec%name
         If(Myid==Master) Call DumpGrid(TempGrid)

        
          Call MPI_Reduce(powerion%Value,TempGrid%Value,Nxy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)   
    
         TempGrid%Value=TempGrid%Value/Dble(Timer) 
         TempGrid%name=powerion%name
         If(Myid==Master) Call DumpGrid(TempGrid)

         ExtotalGrid%Value=ExtotalGrid%Value/Dble(Timer) 
         If(Myid==Master) Call DumpGrid(ExtotalGrid)

       
         Call MPI_Reduce(PowerDensity%Value,TempGrid2%Value,Nxy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)   
        !PowerDensity%Value=PowerDensity%Value+ExGrid%Value
         If(Timer/=0) then
            TempGrid2%Value=(TempGrid2%Value)/Dble(Timer)+ExtotalGrid%Value
         End if
        If(Myid==Master) Call DumpGrid(TempGrid2)  
       !call areapow(TempGrid2%Value,NY,NX,dx)
      
        Return 
     End Subroutine CalPowerDensityInit     

 Subroutine   WeightingPowerDensity(InputParticle,Nx,Ny,Ex,Ey,PowerDensity,dx)
                implicit none
                Type(ParticleBundle),intent(in) :: InputParticle
                Real(8),intent(in) :: dx
                Integer(4),intent(in) :: Nx,Ny
                Real(8),intent(in) :: Ex(Nx,Ny),Ey(Nx,Ny)
                Real(8),intent(out) :: PowerDensity(Nx,Ny)
                Real(8) :: Z1,Z2,Rp1,Rp2,Rn1,Rn2,FactorR,Rho(Nx,Ny)
                Real(8) ::  E1,E2,F11,F12,F21,F22
                Real(8) :: X,Y,CosA,SinA,X1,X2,Y1,Y2
                Real(8) :: Px,Py,Pz,Volume
                Integer(4) :: i,j,N1,N2

                PowerDensity=0.d0
                
                do i=1,InputParticle%Npar
                  
                   
                   N1=Ceiling(InputParticle%PhaseSpace(i)%X)
                   X1=Dble(N1)-InputParticle%PhaseSpace(i)%X
                   X2=1.d0-X1
                   N2=Ceiling(InputParticle%PhaseSpace(i)%Y) 
                   Y1=Dble(N2)-InputParticle%PhaseSpace(i)%Y
                   Y2=1.d0-Y1
                   Rho(N1,N2)=Rho(N1,N2)+X1*Y1
                   Rho(N1,N2+1)=Rho(N1,N2+1)+X1*Y2
                   Rho(N1+1,N2)=Rho(N1+1,N2)+X2*Y1
                   Rho(N1+1,N2+1)=Rho(N1+1,N2+1)+X2*Y2

                   F11=X1*Y1
                   F12=X1*Y2
                   F21=X2*Y1
                   F22=X2*Y2
                        
                   E1=Ex(N1,N2)*F11+Ex(N1+1,N2)*F21+Ex(N1,N2+1)*F12+Ex(N1+1,N2+1)*F22
                   E2=Ey(N1,N2)*F11+Ey(N1+1,N2)*F21+Ey(N1,N2+1)*F12+Ey(N1+1,N2+1)*F22

                   Px=InputParticle%PhaseSpace(i)%Vx*E1
                   Py=InputParticle%PhaseSpace(i)%Vy*E2

                   PowerDensity(N1,N2)=PowerDensity(N1,N2)+F11*(Px+Py)
                   PowerDensity(N1,N2+1)=PowerDensity(N1,N2+1)+F12*(Px+Py)
                   PowerDensity(N1+1,N2)=PowerDensity(N1+1,N2)+F21*(Px+Py)
                   PowerDensity(N1+1,N2+1)=PowerDensity(N1+1,N2+1)+F22*(Px+Py)
           

           end do
          
                 !Volume=1.d0/(dx*dx)
                 PowerDensity(1:Nx,1)=ElectronCharge*PowerDensity(1:Nx,1)
                 do i=2,Ny
                    !  Volume=1.d0/(dx*dx)
                      PowerDensity(1:Nx,i)=ElectronCharge*PowerDensity(1:Nx,i)
                      
                 end do
     return
     end subroutine  WeightingPowerDensity

 Subroutine WeightingEEDFInit()
                implicit none
                Type(ParticleBundle),intent(in) :: InputParticle
                Integer(4),parameter :: Nenergy=1001
                Real(8),parameter :: De=0.1d0
                Type(Grid),save :: EEDF,TempGrid
                Real(8) :: Energy,Eall
                Integer(4),save :: i,j,Timer,EIndex
                Timer=0
                EEDF%Nx=1
                EEDF%Nx=Nenergy
                EEDF%Name='EEPF' 
                TempGrid%Nx=1
                TempGrid%Ny=Nenergy
                TempGrid%Name='EEPF'
                 Return   
         Entry  WeightingEEDF(InputParticle)
                   Timer=Timer+1
                   TempGrid%Value=0.d0
                   !Write(*,*) InputParticle%Mass 
                   do i=1,InputParticle%Npar
                        Call CalEnergy(InputParticle%Mass,InputParticle%PhaseSpace(i),Energy)
                        Energy=Energy*InputParticle%VFactor*InputParticle%VFactor
                        EIndex=Ceiling(Energy/De)
                        If (EIndex<=Nenergy)  then
                             TempGrid% Value(EIndex)=TempGrid% Value(EIndex)+1
                        End if
                    End do
                    EEDF%Value(1:Nenergy)=EEDF%Value(1:Nenergy)+ TempGrid%Value(1:Nenergy)
                   !Write(*,*)  InputParticle%VFactor,Energy
                    Return  
          Entry  WeightingEEDFFinal()
          Call MPI_Reduce(EEDF%Value,TempGrid%Value,Nenergy,mpi_double_precision,MPI_SUM,master,MPI_COMM_WORLD,Ierr)   
          If(Timer/=0) then
               TempGrid%Value=TempGrid%Value/Dble(Timer) 
          End if
do i=1,Nenergy
TempGrid%Value(i)=TempGrid%Value(i)/(i*De)
end do
do i=1,Nenergy
Eall=Eall+TempGrid%Value(i)

end do
TempGrid%Value=TempGrid%Value/Eall
       If(Myid==Master) Call DumpGrid(TempGrid)    
     return
      end subroutine  WeightingEEDFInit

       subroutine DiagDumpParticleInit(InputParticle,InputXMin,InputXMax,InputYmin,InputYmax)
            Implicit none 
            Real(8),intent(in) :: InputXMin,InputXMax,InputYmin,InputYmax
            Real(8),save :: Xmin,Xmax,Ymin,Ymax
            Type(ParticleBundle),intent(inout) :: InputParticle
            Type(ParticleBundle),save :: TempParticle,TempParticle2
            Integer(4):: i 
            Xmin=InputXmin
            Xmax=InputXmax
            Ymin=InputYmin
            Ymax=InputYmax
            Write(TempParticle%Name,*) trim(InputParticle%Name),'X'
            TempParticle%NPar=0
            TempParticle%XFactor=InputParticle%XFactor
            TempParticle%VFactor=InputParticle%VFactor
            TempParticle%Charge=InputParticle%Charge
            TempParticle%Mass=InputParticle%Mass
            Write(TempParticle2%Name,*) trim(InputParticle%Name),'Y' 
            TempParticle2%NPar=0
            TempParticle2%XFactor=InputParticle%XFactor
            TempParticle2%VFactor=InputParticle%VFactor
            TempParticle2%Charge=InputParticle%Charge
            TempParticle2%Mass=InputParticle%Mass
            Return 
       Entry  DiagDumpParticle(InputParticle)
             do i=InputParticle%NPar,1,-1

                    
                            If(InputParticle%PhaseSpace(i)%X< XMIN) then
                                           Call AddParticle(InputParticle%PhaseSpace(i),TempParticle)
              
                             Else if(InputParticle%PhaseSpace(i)%Y<Ymin)  then
                                           Call AddParticle(InputParticle%PhaseSpace(i),TempParticle2)
                              
                             End if
             End do
 
             Return 
       Entry  DiagDumpParticleFinal()
             Write(*,*)  TempParticle%NPar,TempParticle2%NPar,Myid
             Call ParallelDumpParticle(TempParticle)
             Call ParallelDumpParticle(TempParticle2)
             Return
       End subroutine DiagDumpParticleInit
End Module Diagnostics