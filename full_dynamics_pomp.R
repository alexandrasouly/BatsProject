# full dynamics SIR/SIRS/SILI
# deterministic equations
# same equs as in Aaron's supp material
# maternal immunity included
# here for now E and R gives maternal immunity as in supp eqs, not github equs

library(pomp)

det_model <- Csnippet("

Sbirths =b*(Sf +If);
MaBirths =b*(Rf + Ef);

betSN =omega_m+mj*(N/kappa) + beta*(If+Im+Ij+In);
betSj=mu+mj*(N/kappa)+(beta)*(If+Im+In+Ij);
betSmf=m*(N/kappa)+(beta)*(If+Im+In+Ij);
betI=beta*(If+Im+In+Ij);

#newborn transitions

if (Ma+ (MaBirths) - (omega_m+mj*(N/kappa))*Ma <=0) 
{DMa =0;}
else
{DMa=Ma+ (MaBirths) - (omega_m+mj*(N/kappa))*Ma;}


if (Sn + (Sbirths) - (betSN*Sn) + omega*Rn <=0)
  {DSn = 0;}
else
  {DSn = Sn + (Sbirths) - (betSN*Sn) + omega*Rn;}
  
  
if (En - ((omega_m+mj*(N/kappa)+epsilon)*En)+rho*In <=0)
  {DEn =0;}
else 
  {DEn = En - ((omega_m+mj*(N/kappa)+epsilon)*En)+rho*In;}
  

if (In -((omega_m+mj*(N/kappa)+gamma+rho)*In)+(betI*Sn)+epsilon*En <=0)
 {DIn = 0;}
else 
  {DIn = In -((omega_m+mj*(N/kappa)+gamma+rho)*In)+(betI*Sn)+epsilon*En;}

if (Rn -(omega_m+mj*(N/kappa)+omega)*Rn+gamma*In <=0)
  {DRn=0;}
else
  {DRn = Rn -(omega_m+mj*(N/kappa)+omega)*Rn +gamma*In;}                      
 
 
## transitions between juvenile compartments:

if (Sj+omega_m*(Sn+Ma) -betSj*Sj +omega*Rj <=0) 
  {DSj =0;}
else 
  {DSj =Sj+omega_m*(Sn+Ma) -betSj*Sj+omega*Rj;}

if (Ej+omega_m*En-(mu+mj*(N/kappa)+epsilon)*Ej+rho*Ij <=0) 
{DEj =0;} 
else 
{DEj =  Ej+omega_m*En-(mu+mj*(N/kappa)+epsilon)*Ej+rho*Ij;}

if (Ij+omega_m*In-(mu+mj*(N/kappa)+gamma+rho)*Ij+betI*Sj+epsilon*Ej <=0) 
{DIj =0;}
else 
{DIj =Ij+omega_m*In-(mu+mj*(N/kappa)+gamma+rho)*Ij+betI*Sj+epsilon*Ej;}

if (Rj+omega_m*Rn-(mu+mj*(N/kappa)+omega)*Rj+gamma*Ij <=0) 
{DRj =0;}
else 
{DRj = Rj+omega_m*Rn-(mu+mj*(N/kappa)+omega)*Rj+gamma*Ij;}

## transitions between adult male compartments:

if (Sm+ mu*(Sj/2)-betSmf*Sm+omega*Rm <=0)
{DSm =0;}
else 
{DSm =Sm+ mu*(Sj/2)-betSmf*Sm+omega*Rm;}

if (Em+ mu*(Ej/2) -(m*(N/kappa)+epsilon)*Em+rho*Im <=0) 
{DEm =0;}
else
{DEm =Em+ mu*(Ej/2) -(m*(N/kappa)+epsilon)*Em+rho*Im;}


if (Im+ mu*(Ij/2)-(m*(N/kappa)+gamma+rho)*Im+(betI*Sm)+epsilon*Em <=0)
{DIm =0;}
else 
{DIm =Im+ mu*(Ij/2)-(m*(N/kappa)+gamma+rho)*Im+(betI*Sm)+epsilon*Em;}

if (Rm+ mu*(Rj/2)-(m*(N/kappa)+omega)*Rm+gamma*Im <=0)
{DRm =0;}
else
{DRm = Rm+ mu*(Rj/2)-(m*(N/kappa)+omega)*Rm+gamma*Im;}

## transitions between adult female compartments:
if (Sf+ mu*(Sj/2)-betSmf*Sf+omega*Rf <=0)
{DSf =0;}
else 
{DSf =Sf+ mu*(Sj/2)-betSmf*Sf+omega*Rf;}

if (Ef+ mu*(Ej/2) -(m*(N/kappa)+epsilon)*Ef+rho*If <=0) 
{DEf =0;}
else
{DEf =Ef+ mu*(Ej/2) -(m*(N/kappa)+epsilon)*Ef+rho*If;}


if (If+ mu*(Ij/2)-(m*(N/kappa)+gamma+rho)*If+(betI*Sf)+epsilon*Ef <=0)
{DIf =0;}
else 
{DIf =If+ mu*(Ij/2)-(m*(N/kappa)+gamma+rho)*If+(betI*Sf)+epsilon*Ef;}

if (Rf+ mu*(Rj/2)-(m*(N/kappa)+omega)*Rf+gamma*If <=0)
{DRf =0;}
else
{DRf = Rf+ mu*(Rj/2)-(m*(N/kappa)+omega)*Rf+gamma*If;}


")