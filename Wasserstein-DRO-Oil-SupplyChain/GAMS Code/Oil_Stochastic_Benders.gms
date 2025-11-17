sets
i /1*3/
j /1*5/
k /1*3/
m /1*15/
t /1*12/
sa /1*10/
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
ca4(i,k,sa)
ca4av(i,k)
ca4_upper(i,k)
ca4_lower(i,k)
d(i,m,t,sa)
dav(i,m,t)
d_upper(i,m,t)
d_lower(i,m,t)
alpha(i)
u(j,jp)
BigM /1000/
theta /0/
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

*/////////////////////////Dual Sub Problem /////////////////////////////////////
parameters
x_star(i,j,k)
zeta_star(i,j,n);
x_star(i,j,k) = 10;
zeta_star(i,j,n) = 0;
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
pivar(i,k,t);
Positive Variables eta, nu, rho, pivar;

Variables
zdsp
gammavar(i,k,t)
mu(i,k,t)
eta(i,k,t)
nu(i,k,t)
rho(i,m,t)
pivar(i,k,t);
Free Variables zdsp, gammavar, mu;
Positive Variables eta, nu, rho, pivar;

Equations
objdsp
cons1dsp(i,k,t)
cons2dsp(i,k,m,t)
cons3dsp(i,k,t)
cons4dsp(i,k,t)
cons5dsp(i,k,t)
cons6dsp(i,k,t);

objdsp .. zdsp =e= (-1)*sum((i,k,t),ca4_hat(i,k)*nu(i,k,t)) - sum((i,k,t)$(ord(t)=1),(ca2(i,k)-sum(j,x_star(i,j,k)))*eta(i,k,t)) - sum((i,k,t)$(ord(t)>1),ca2(i,k)*eta(i,k,t))
                    - sum((i,k,j,t)$(ord(t)=1),x_star(i,j,k)*gammavar(i,k,t)) - sum((i,k,t)$(ord(t)<card(t)),ca3(i,k)*pivar(i,k,t))
                    - sum((i,m,t),d_hat(i,m,t)*rho(i,m,t));

cons1dsp(i,k,t) .. (-1)*alpha(i)*gammavar(i,k,t) + mu(i,k,t) - nu(i,k,t) =l= c2(i);
cons2dsp(i,k,m,t) .. (-1)*mu(i,k,t) - rho(i,m,t) =l= c3(i,k,m) - v1(i,t);
cons3dsp(i,k,t) .. gammavar(i,k,t) + eta(i,k,t) =l= c4(i,t);
cons4dsp(i,k,t) .. gammavar(i,k,t) =g= v2(i,t);
cons5dsp(i,k,t)$(ord(t)<card(t)) .. (-1)*gammavar(i,k,t) + gammavar(i,k,t+1) + eta(i,k,t+1) =l= c5(i,k,t);
cons6dsp(i,k,t)$(ord(t)<card(t)) .. (-1)*mu(i,k,t) + mu(i,k,t+1) + pivar(i,k,t) =l= c6(i,k,t);

Model dualsubproblem
/objdsp,
cons1dsp,
cons2dsp,
cons3dsp,
cons4dsp,
cons5dsp,
cons6dsp/;

*/////////////////////////Master Problem /////////////////////////////////////
Set optcut(sa,iter);
optcut(sa,iter) = no;
Variables
lb
phi_cut(sa)
zeta(i,j,n)
x(i,j,k)
delta(j);
Sos2 Variables zeta;
Free Variables lb, phi_cut(sa);
Positive Variables x;
Binary Variables delta;
Parameters
gammavar_star(i,k,t,sa,iter)
mu_star(i,k,t,sa,iter)
eta_star(i,k,t,sa,iter)
nu_star(i,k,t,sa,iter)
rho_star(i,m,t,sa,iter)
pivar_star(i,k,t,sa,iter);

Equations
objmaster
cons1master(i,j)
cons2master(i,k)
cons3master(i,j)
cons4master(j,jp)
cons5master(i,j,k)
cons6master(i,j)
optimalitycut(sa,iter);

objmaster .. lb =e= sum((i,j,n),f1(i,j,n)*zeta(i,j,n)) + sum((i,j,k),c1(i,j,k)*x(i,j,k))  + sampel_size*sum(sa,phi_cut(sa));
cons1master(i,j) .. sum(k,x(i,j,k)) =l= ca1(i,j);
cons2master(i,k) .. sum(j,x(i,j,k)) =l= ca2(i,k);
cons3master(i,j) .. sum(k,x(i,j,k)) =e= sum(n,g(i,j,n)*zeta(i,j,n));
cons4master(j,jp) .. delta(j) + delta(jp) =l= 1 + u(j,jp);
cons5master(i,j,k) .. x(i,j,k) =l= BigM*delta(j);
cons6master(i,j) .. sum(n,zeta(i,j,n)) =e= 1;
optimalitycut(sa,iter1)$(optcut(sa,iter1)) .. phi_cut(sa) =g= (-1)*sum((i,k,t),ca4_hat(i,k)*nu_star(i,k,t,sa,iter1)) - sum((i,k,t)$(ord(t)=1),(ca2(i,k)-sum(j,x(i,j,k)))*eta_star(i,k,t,sa,iter1)) - sum((i,k,t)$(ord(t)>1),ca2(i,k)*eta_star(i,k,t,sa,iter1))
                    - sum((i,k,j,t)$(ord(t)=1),x(i,j,k)*gammavar_star(i,k,t,sa,iter1)) - sum((i,k,t)$(ord(t)<card(t)),ca3(i,k)*pivar_star(i,k,t,sa,iter1))
                    - sum((i,m,t),d_hat(i,m,t)*rho_star(i,m,t,sa,iter1));


Model masterproblem /objmaster, cons1master, cons2master, cons3master, cons4master,
                     cons5master, cons6master, optimalitycut/;

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
loop(sa,
ca4_hat(i,k) = ca4(i,k,sa);
d_hat(i,m,t) = d(i,m,t,sa);
option optcr=0;
option mip=cplex;
solve dualsubproblem maximizing zdsp using mip;
ub(iter)=ub(iter)+zdsp.l;
optcut(sa,iter)=yes;
gammavar_star(i,k,t,sa,iter) = gammavar.l(i,k,t);
mu_star(i,k,t,sa,iter) = mu.l(i,k,t);
eta_star(i,k,t,sa,iter) = eta.l(i,k,t);
nu_star(i,k,t,sa,iter) = nu.l(i,k,t);
rho_star(i,m,t,sa,iter) = rho.l(i,m,t);
pivar_star(i,k,t,sa,iter) = pivar.l(i,k,t);
);
uperbound=sum((i,j,n),f1(i,j,n)*zeta_star(i,j,n)) + sum((i,j,k),c1(i,j,k)*x_star(i,j,k)) + sampel_size*ub(iter)
option optcr=0;
option mip=cplex;
solve masterproblem minimizing lb using mip;
x_star(i,j,k) = x.l(i,j,k);
zeta_star(i,j,n) = zeta.l(i,j,n);
display "Lower Bound Value is: ------------------------------------------------ ",lb.l;
abort$(masterproblem.modelstat<>1) "original model not solved optimality";
lowerbound=lb.l;
result(iter,'lb')=lowerbound;
result(iter,'ub')=uperbound;
result(iter,'phicut') = sum(sa,phi_cut.l(sa));
nonconverged$((uperbound-lowerbound)<0.1)=no;
);
endTime = jnow;
executionTime = (startTime-endTime)*86400;
display uperbound, lowerbound, result, x.l, zeta.l, delta.l, executionTime, x_star, ub;
