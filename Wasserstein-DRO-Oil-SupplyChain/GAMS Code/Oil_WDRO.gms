sets
i /1*3/
j /1*5/
k /1*3/
m /1*15/
t /1*12/
s /1*10/
n /1*61/
iter /1*5/;
Alias (iter,iter1);
Alias (j, jp);
parameters
f1(i,j,n)
g(i,j,n)
c1(i,j,k)
c2(i)
c3(i,k,m)
c4(i,t)
c5(i,k,t)
c6(i,k,t)
v1(i,t)
v2(i,t)
ca1(i,j)
ca2(i,k)
ca3(i,k)
ca4(i,k,s)
ca4_upper(i,k)
ca4_lower(i,k)
d(i,m,t,s)
d_upper(i,m,t)
d_lower(i,m,t)
alpha(i)
u(j,jp)
BigM /1000000/
theta /10000/
sampel_size /0.1/;
*////////////////////////// Include Parameter Data ///////////////////////
$call gdxxrw oildata.xlsx  par=f1 rng=f1! rdim=3 cdim=0
$gdxin oildata.gdx
$load f1
$gdxin

$call gdxxrw oildata.xlsx  par=g rng=g! rdim=3 cdim=0
$gdxin oildata.gdx
$load g
$gdxin

$call gdxxrw oildata.xlsx  par=c1 rng=c1! rdim=3 cdim=0
$gdxin oildata.gdx
$load c1
$gdxin

$call gdxxrw oildata.xlsx  par=c2 rng=c2! rdim=1 cdim=0
$gdxin oildata.gdx
$load c2
$gdxin

$call gdxxrw oildata.xlsx  par=c3 rng=c3! rdim=3 cdim=0
$gdxin oildata.gdx
$load c3
$gdxIn

$call gdxxrw oildata.xlsx  par=c4 rng=c4! rdim=2 cdim=0
$gdxin oildata.gdx
$load c4
$gdxIn

$call gdxxrw oildata.xlsx  par=c5 rng=c5! rdim=3 cdim=0
$gdxin oildata.gdx
$load c5
$gdxIn

$call gdxxrw oildata.xlsx  par=c6 rng=c6! rdim=3 cdim=0
$gdxin oildata.gdx
$load c6
$gdxIn

$call gdxxrw oildata.xlsx  par=v1 rng=v1! rdim=2 cdim=0
$gdxin oildata.gdx
$load v1
$gdxIn

$call gdxxrw oildata.xlsx  par=v2 rng=v2! rdim=2 cdim=0
$gdxin oildata.gdx
$load v2
$gdxIn

$call gdxxrw oildata.xlsx  par=ca1 rng=ca1! rdim=2 cdim=0
$gdxin oildata.gdx
$load ca1
$gdxIn

$call gdxxrw oildata.xlsx  par=ca2 rng=ca2! rdim=2 cdim=0
$gdxin oildata.gdx
$load ca2
$gdxIn

$call gdxxrw oildata.xlsx  par=ca3 rng=ca3! rdim=2 cdim=0
$gdxin oildata.gdx
$load ca3
$gdxIn

$call gdxxrw oildata.xlsx  par=ca4 rng=ca4! rdim=3 cdim=0
$gdxin oildata.gdx
$load ca4
$gdxIn

$call gdxxrw oildata.xlsx  par=ca4_upper rng=ca4_upper! rdim=2 cdim=0
$gdxin oildata.gdx
$load ca4_upper
$gdxIn

$call gdxxrw oildata.xlsx  par=ca4_lower rng=ca4_lower! rdim=2 cdim=0
$gdxin oildata.gdx
$load ca4_lower
$gdxIn

$call gdxxrw oildata.xlsx  par=d rng=d! rdim=4 cdim=0
$gdxin oildata.gdx
$load d
$gdxIn

$call gdxxrw oildata.xlsx  par=d_upper rng=d_upper! rdim=3 cdim=0
$gdxin oildata.gdx
$load d_upper
$gdxIn

$call gdxxrw oildata.xlsx  par=d_lower rng=d_lower! rdim=3 cdim=0
$gdxin oildata.gdx
$load d_lower
$gdxIn

$call gdxxrw oildata.xlsx  par=alpha rng=alpha! rdim=1 cdim=0
$gdxin oildata.gdx
$load alpha
$gdxIn

$call gdxxrw oildata.xlsx  par=u rng=u! rdim=2 cdim=0
$gdxin oildata.gdx
$load u
$gdxIn
*///////////////////////// First Stage Variabels ///////////////////////////////
variables
x(i,j,k)
zeta(i,j,n)
delta(j)
landa
phi(s);
Binary Variables delta;
Positive Variable x, zeta, landa;
Free Variable phi;

*/////////////////////////Dual Sub Problem /////////////////////////////////////
parameters
x_star(i,j,k)
zeta_star(i,j,n)
landa_star;
x_star(i,j,k) = 0;
zeta_star(i,j,n) = 0;
landa_star = 100;
parameters
ca4_hat(i,k)
d_hat(i,m,t);
Variables
zdsp
gammavar(i,k,t)
mu(i,k,t)
eta(i,k,t)
nu(i,k,t)
rho(i,m,t)
pivar(i,k,t)
omega_plus(i,k,t)
omega_minus(i,k,t)
comega_plus(i,m,t)
comega_minus(i,m,t)
chi_plus(i,k,t)
chi_minus(i,k,t)
gamma_plus(i,m,t)
gamma_minus(i,m,t);
Positive Variables eta, nu, rho, pivar, omega_plus, omega_minus, comega_plus, comega_minus;
Binary Variables chi_plus, chi_minus, gamma_plus, gamma_minus;

Equations
objdsp
cons1dsp(i,k,t)
cons2dsp(i,k,m,t)
cons3dsp(i,k,t)
cons4dsp(i,k,t)
cons5dsp(i,k,t)
cons6dsp(i,k,t)
cons7dsp(i,k,t)
cons8dsp(i,k,t)
cons9dsp(i,k,t)
cons10dsp(i,k,t)
cons11dsp(i,k,t)
cons12dsp(i,k,t)
cons13dsp(i,k,t)
cons14dsp(i,m,t)
cons15dsp(i,m,t)
cons16dsp(i,m,t)
cons17dsp(i,m,t)
cons18dsp(i,m,t)
cons19dsp(i,m,t)
cons20dsp(i,m,t);

objdsp .. zdsp =e= (-1)*sum((i,k,t),(ca4_upper(i,k)-ca4_hat(i,k))*omega_plus(i,k,t) + (ca4_lower(i,k)-ca4_hat(i,k))*omega_minus(i,k,t) + ca4_hat(i,k)*nu(i,k,t))
                    - sum((i,k,t)$(ord(t)=1),(ca2(i,k)-sum(j,x_star(i,j,k)))*eta(i,k,t)) - sum((i,k,t)$(ord(t)>1),ca2(i,k)*eta(i,k,t))
                    - sum((i,k,j,t)$(ord(t)=1),x_star(i,j,k)*gammavar(i,k,t)) - sum((i,k,t)$(ord(t)<card(t)),ca3(i,k)*pivar(i,k,t))
                    - sum((i,m,t),(d_upper(i,m,t)-d_hat(i,m,t))*comega_plus(i,m,t) + (d_lower(i,m,t)-d_hat(i,m,t))*comega_minus(i,m,t) + d_hat(i,m,t)*rho(i,m,t))
                    - landa_star*sum((i,k,t),(ca4_upper(i,k)-ca4_hat(i,k))*chi_plus(i,k,t) + (ca4_hat(i,k)-ca4_lower(i,k))*chi_minus(i,k,t))
                    - landa_star*sum((i,m,t),(d_upper(i,m,t)-d_hat(i,m,t))*gamma_plus(i,m,t) + (d_hat(i,m,t)-d_lower(i,m,t))*gamma_minus(i,m,t));

cons1dsp(i,k,t) .. (-1)*alpha(i)*gammavar(i,k,t) + mu(i,k,t) - nu(i,k,t) =l= c2(i);
cons2dsp(i,k,m,t) .. (-1)*mu(i,k,t) - rho(i,m,t) =l= c3(i,k,m) - v1(i,t);
cons3dsp(i,k,t) .. gammavar(i,k,t) + eta(i,k,t) =l= c4(i,t);
cons4dsp(i,k,t) .. gammavar(i,k,t) =g= v2(i,t);
cons5dsp(i,k,t)$(ord(t)<card(t)) .. (-1)*gammavar(i,k,t) + gammavar(i,k,t+1) + eta(i,k,t+1) =l= c5(i,k,t);
cons6dsp(i,k,t)$(ord(t)<card(t)) .. (-1)*mu(i,k,t) + mu(i,k,t+1) + pivar(i,k,t) =l= c6(i,k,t);
cons7dsp(i,k,t) .. omega_plus(i,k,t) =l= BigM*chi_plus(i,k,t);
cons8dsp(i,k,t) .. omega_plus(i,k,t) =l= nu(i,k,t);
cons9dsp(i,k,t) .. omega_plus(i,k,t) =g= nu(i,k,t) - BigM*(1-chi_plus(i,k,t));
cons10dsp(i,k,t) .. omega_minus(i,k,t) =l= BigM*chi_minus(i,k,t);
cons11dsp(i,k,t) .. omega_minus(i,k,t) =l= nu(i,k,t);
cons12dsp(i,k,t) .. omega_minus(i,k,t) =g= nu(i,k,t) - BigM*(1-chi_minus(i,k,t));
cons13dsp(i,k,t) .. chi_plus(i,k,t) + chi_minus(i,k,t) =l= 1;
cons14dsp(i,m,t) .. comega_plus(i,m,t) =l= BigM*gamma_plus(i,m,t);
cons15dsp(i,m,t) .. comega_plus(i,m,t) =l= rho(i,m,t);
cons16dsp(i,m,t) .. comega_plus(i,m,t) =g= rho(i,m,t) - BigM*(1-gamma_plus(i,m,t));
cons17dsp(i,m,t) .. comega_minus(i,m,t) =l= BigM*gamma_minus(i,m,t);
cons18dsp(i,m,t) .. comega_minus(i,m,t) =l= rho(i,m,t);
cons19dsp(i,m,t) .. comega_minus(i,m,t) =g= rho(i,m,t) - BigM*(1-gamma_minus(i,m,t));
cons20dsp(i,m,t) .. gamma_plus(i,m,t) + gamma_minus(i,m,t) =l= 1;

Model dualsubproblem
/objdsp,
cons1dsp,
cons2dsp,
cons3dsp,
cons4dsp,
cons5dsp,
cons6dsp,
cons7dsp,
cons8dsp,
cons9dsp,
cons10dsp,
cons11dsp,
cons12dsp,
cons13dsp,
cons14dsp,
cons15dsp,
cons16dsp,
cons17dsp,
cons18dsp,
cons19dsp,
cons20dsp/;
*///////////////////////// Master Problem //////////////////////////////////////
Set optcut(s,iter);
optcut(s,iter) = no;
Variables
phi_cut(s)
lb;
Positive Variables phi_cut;
Parameters
gammavar_star(i,k,t,s,iter)
mu_star(i,k,t,s,iter)
eta_star(i,k,t,s,iter)
nu_star(i,k,t,s,iter)
rho_star(i,m,t,s,iter)
pivar_star(i,k,t,s,iter)
omega_plus_star(i,k,t,s,iter)
omega_minus_star(i,k,t,s,iter)
comega_plus_star(i,m,t,s,iter)
comega_minus_star(i,m,t,s,iter)
chi_plus_star(i,k,t,s,iter)
chi_minus_star(i,k,t,s,iter)
gamma_plus_star(i,m,t,s,iter)
gamma_minus_star(i,m,t,s,iter);
Equations
objmaster
cons1master(i,j)
cons2master(i,k)
cons3master(i,j)
cons4master(j,jp)
cons5master(i,j,k)
cons6master(i,j)
optimalitycut(s,iter);

objmaster .. lb =e= sum((i,j,n),f1(i,j,n)*zeta(i,j,n)) + sum((i,j,k),c1(i,j,k)*x(i,j,k)) + landa*theta + sampel_size*sum(s,phi_cut(s));
cons1master(i,j) .. sum(k,x(i,j,k)) =l= ca1(i,j);
cons2master(i,k) .. sum(j,x(i,j,k)) =l= ca2(i,k);
cons3master(i,j) .. sum(k,x(i,j,k)) =e= sum(n,g(i,j,n)*zeta(i,j,n));
cons4master(j,jp) .. delta(j) + delta(jp) =l= 1 + u(j,jp);
cons5master(i,j,k) .. x(i,j,k) =l= BigM*delta(j);
cons6master(i,j) .. sum(n,zeta(i,j,n)) =e= 1;
optimalitycut(s,iter1)$(optcut(s,iter1)) .. phi_cut(s) =g= (-1)*sum((i,k,t),(ca4_upper(i,k)-ca4_hat(i,k))*omega_plus_star(i,k,t,s,iter1) + (ca4_lower(i,k)-ca4_hat(i,k))*omega_minus_star(i,k,t,s,iter1) + ca4_hat(i,k)*nu_star(i,k,t,s,iter1))
                    - sum((i,k,t)$(ord(t)=1),(ca2(i,k)-sum(j,x(i,j,k)))*eta_star(i,k,t,s,iter1)) - sum((i,k,t)$(ord(t)>1),ca2(i,k)*eta_star(i,k,t,s,iter1))
                    - sum((i,k,j,t)$(ord(t)=1),x(i,j,k)*gammavar_star(i,k,t,s,iter1)) - sum((i,k,t)$(ord(t)<card(t)),ca3(i,k)*pivar_star(i,k,t,s,iter1))
                    - sum((i,m,t),(d_upper(i,m,t)-d_hat(i,m,t))*comega_plus_star(i,m,t,s,iter1) + (d_lower(i,m,t)-d_hat(i,m,t))*comega_minus_star(i,m,t,s,iter1) + d_hat(i,m,t)*rho_star(i,m,t,s,iter1))
                    - landa*sum((i,k,t),(ca4_upper(i,k)-ca4_hat(i,k))*chi_plus_star(i,k,t,s,iter1) + (ca4_hat(i,k)-ca4_lower(i,k))*chi_minus_star(i,k,t,s,iter1))
                    - landa*sum((i,m,t),(d_upper(i,m,t)-d_hat(i,m,t))*gamma_plus_star(i,m,t,s,iter1) + (d_hat(i,m,t)-d_lower(i,m,t))*gamma_minus_star(i,m,t,s,iter1));

Model master /objmaster,cons1master,cons2master,cons3master,cons4master,cons5master,cons6master,optimalitycut/;
*////////////////////////// Benders Algorithm //////////////////////////////////
Parameter result(iter,*);
Parameter nonconverged /yes/;
Parameter ub(iter);
ub(iter)=0;
Scalar uperbound /+inf/;
Scalar lowerbound /-inf/;
Scalars startTime, endTime, executionTime;

startTime = jnow;
loop(iter$(nonconverged),
loop(s,
ca4_hat(i,k) = ca4(i,k,s);
d_hat(i,m,t) = d(i,m,t,s);
option optcr=0;
option mip=cplex;
solve dualsubproblem maximizing zdsp using mip;
ub(iter)=ub(iter)+zdsp.l;
gammavar_star(i,k,t,s,iter) = gammavar.l(i,k,t);
mu_star(i,k,t,s,iter) = mu.l(i,k,t);
eta_star(i,k,t,s,iter) = eta.l(i,k,t);
nu_star(i,k,t,s,iter) = nu.l(i,k,t);
rho_star(i,m,t,s,iter) = rho.l(i,m,t);
pivar_star(i,k,t,s,iter) = pivar.l(i,k,t);
omega_plus_star(i,k,t,s,iter) = omega_plus.l(i,k,t);
omega_minus_star(i,k,t,s,iter) = omega_minus.l(i,k,t);
comega_plus_star(i,m,t,s,iter) = comega_plus.l(i,m,t);
comega_minus_star(i,m,t,s,iter) = comega_minus.l(i,m,t);
chi_plus_star(i,k,t,s,iter) = chi_plus.l(i,k,t);
chi_minus_star(i,k,t,s,iter) = chi_minus.l(i,k,t);
gamma_plus_star(i,m,t,s,iter) = gamma_plus.l(i,m,t);
gamma_minus_star(i,m,t,s,iter) = gamma_minus.l(i,m,t);
optcut(s,iter)=yes;
);
uperbound=sum((i,j,n),f1(i,j,n)*zeta_star(i,j,n)) + sum((i,j,k),c1(i,j,k)*x_star(i,j,k)) + landa_star*theta + sampel_size*ub(iter)
option optcr=0;
option mip=cplex;
solve master minimizing lb using mip;
x_star(i,j,k) = x.l(i,j,k);
zeta_star(i,j,n) = zeta.l(i,j,n);
landa_star = landa.l;
display "Lower Bound Value is: ------------------------------------------------ ",lb.l;
abort$(master.modelstat<>1) "original model not solved optimality";
lowerbound=lb.l;
result(iter,'lb')=lowerbound;
result(iter,'ub')=uperbound;
result(iter,'cut_lhs')=sum(s,phi_cut.l(s));
nonconverged$((uperbound-lowerbound)<0.1)=no;
);
endTime = jnow;
executionTime = (startTime-endTime)*86400;
display uperbound, lowerbound, result, x.l, zeta.l, delta.l, landa.l, executionTime, optcut;
