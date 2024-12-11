Module MaxwellModule
    Use TypeModule
    USE Constants
    Use EnergyKai
    Implicit none
	contains
   
subroutine Maxwellian2(InputGas,InputParticle)
       implicit none
       Type(Gas),intent(in) :: InputGas
       Type(ParticleOne),intent(out) :: InputParticle
       real(8) ::  Mass,Temprature 
       real(8) :: V,Beta,FuncA,FuncB
       real(8) :: Theta,CosTheta,SinTheta,Fai 
       Mass=InputGas%MGas
       Temprature=InputGas%TGas
       Beta=1.d0/(kB*Temprature)
       FuncA=1.d0
       FuncB=0.d0
      do while(FuncA>FuncB)
         Call DRandom(R)
              FuncA=R*R
         Call DRandom(R)
              FuncB=-exp*R*Dlog(R)
      end do

      V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)

      Call RandomVelocity(V,InputParticle%Vx,InputParticle%Vy,InputParticle%Vz)
     return 
   end subroutine Maxwellian2

    subroutine Maxwellian(InputGas,InputParticle)
       implicit none
       Type(Gas),intent(in) :: InputGas
       Type(ParticleOne),intent(out) :: InputParticle
       real(8) ::  Mass,Temprature 
       real(8) :: V,Beta,FuncA,FuncB
       real(8) :: Theta,CosTheta,SinTheta,Fai 
       Mass=InputGas%MGas
       Temprature=InputGas%TGas
       Beta=1.d0/(kB*Temprature)
       FuncA=1.d0
       FuncB=0.d0
      do while(FuncA>FuncB)
         Call DRandom(R)
              FuncA=R*R
         Call DRandom(R)
              FuncB=-exp*R*Dlog(R)
      end do
      V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)
      Call RandomVelocity(V,InputParticle%Vx,InputParticle%Vy,InputParticle%Vz)
     return 
   end subroutine Maxwellian
  
   subroutine RandomVelocity(V,Vx,Vy,Vz)
       implicit none
       Real(8),intent(in) ::  V
       Real(8),intent(out) ::  Vx,Vy,Vz  
       Real(8) :: Fai,CosFai,SinFai
       Real(8) :: Theta,CosTheta,FcosTheta,SinTheta
       Call DRandom(R)
       CosTheta=IsotropicCosKai(Theta)
        SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)
        Call DRandom(R)
        Fai=2.d0*PI*R
        Vx=V*CosTheta
        Vy=V*SinTheta*DCos(Fai)
        Vz=V*SinTheta*Dsin(Fai)
       return
   end subroutine RandomVelocity   
  
End Module MaxwellModule   

