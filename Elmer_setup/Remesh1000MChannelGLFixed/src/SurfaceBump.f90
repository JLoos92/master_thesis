FUNCTION SurfaceBump( Model, nodenumber, Time) RESULT(BedBump)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: BedSol
   INTEGER, POINTER :: BedPerm(:)
   REAL(kind=dp), POINTER :: BedVal(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: BedBump, Time, x,y,y0,x0,sigmax,sigmay,Ampl,Bump
   REAL(KIND=dp) :: MaxGrowTime
   LOGICAL :: FirstTime=.True., UnFoundFatal

   SAVE FirstTime
   MaxGrowTime = 10.0
   Ampl = 150.0
   sigmax = 500
   sigmay = 500
   x0 = 1056000
   y0 = 0
   x = Model % Nodes % x(nodenumber)
   y = Model % Nodes % y(nodenumber)
   
   !Get Bedrock without catching any error messages if fields don't exist
   BedSol => VariableGet( Model % Variables, 'Zs',UnFoundFatal=UnFoundFatal)
   BedPerm => BedSol % Perm
   BedVal => BedSol % Values

   !Bump = exp(-((x-x0)**2/(2*sigmax**2)))
   if (Time <= MaxGrowTime) then
      Bump =(Time/MaxGrowTime)*Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   else
      Bump =Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   end if
   if (x>x0) then
      Bump = 0.0
   end if
   !write(*,*) 'BUMP',Bump
   BedBump = BedVal(BedPerm(nodenumber))+Bump

   !Get Distance without catching any error messages if fields don't exist
   !DistanceSol => VariableGet( Model % Variables, 'Distance',UnFoundFatal=UnFoundFatal)
   !DistancePerm => DistanceSol % Perm
   !DistanceVal => DistanceSol % Values
!

  !alpha = 0.5
  !rho = 1-np.exp(-0.0001*distance) #transition GL to Ambient
  !G = 0.001 ## melting away from GL relative to H^alpha
  !A = 0.1   ##melting near GL relative to H^alpha
  !bmb = thickness**(alpha)*(rho*G+(1-rho)*A)


   !alpha = 0.4
   !alpha = 0.76
   !G = 0.00001
   !A = 0.06
   !rho = 1 - exp(-0.0002* DistanceVal(DistancePerm(nodenumber)))

   !BMBOut = BMBMultiplier*DepthVal(DepthPerm(nodenumber))**alpha*(rho*G+(1-rho)*A)
   !print *, BMBOut
   !write(*,*) 'BMBOut:',Model % Nodes % x(nodenumber),Model % Nodes % y(nodenumber),BMBOut

END FUNCTION SurfaceBump
