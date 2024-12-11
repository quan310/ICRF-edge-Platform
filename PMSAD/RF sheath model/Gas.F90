Module GasModule
     Use TypeModule
     Use Constants
     Implicit none
     contains  
              subroutine  GasInit(GasPressure,GasTemperature,InputGasPhysical,OutputGas,Shift)
                  Implicit none
                  Integer(4), intent(inout) :: Shift
                  Real(8), intent(in) :: GasPressure,GasTemperature
                  Type(GasPhysical), intent(in) :: InputGasPhysical 
                  Type(Gas), intent(out) ::  OutputGas
                  OutputGas%Shift=Shift
                  OutputGas%NSpecy=InputGasPhysical%NSpecy
                  OutputGas%MGas=InputGasPhysical%MGas
                  OutputGas%TGas=GasTemperature
                  OutputGas%PGas=GasPressure*mTorrtoPa
                  OutputGas%BetaMax=InputGasPhysical%BetaMax
                  return
             End  subroutine  GasInit
 End Module  GasModule

Module EnergyKai
     Use Constants
     Implicit none
     contains
	 
     Function IsotropicCosKai(Energy) 
       Implicit none
       Real(8) :: IsotropicCosKai,Energy
       Call DRandom(R)
       IsotropicCosKai=1.d0-2.d0*R
       Return
    end Function IsotropicCosKai
             
End Module  EnergyKai