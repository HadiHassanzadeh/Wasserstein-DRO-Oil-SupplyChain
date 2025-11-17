sets
i /1*3/
j /1*5/
k /1*3/
m /1*15/
t /1*12/
s /1*40/
n /1*61/
iter /1*200/;
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
ca4av(i,k)
ca4_upper(i,k)
ca4_lower(i,k)
d(i,m,t,s)
dav(i,m,t)
d_upper(i,m,t)
d_lower(i,m,t)
alpha(i)
u(j,jp)
BigM /1000/
theta /0/
sampel_size /0.025/;
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


ca4av(i,k) = sum(s,ca4(i,k,s))/40;
dav(i,m,t) = sum(s,d(i,m,t,s))/40;
parameters
x_star(i,j,k);


x_star(i,j,k) = 10;

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

objdsp .. zdsp =e= (-1)*sum((i,k,t),ca4av(i,k)*nu(i,k,t)) - sum((i,k,t)$(ord(t)=1),(ca2(i,k)-sum(j,x_star(i,j,k)))*eta(i,k,t)) - sum((i,k,t)$(ord(t)>1),ca2(i,k)*eta(i,k,t))
                    - sum((i,k,j,t)$(ord(t)=1),x_star(i,j,k)*gammavar(i,k,t)) - sum((i,k,t)$(ord(t)<card(t)),ca3(i,k)*pivar(i,k,t))
                    - sum((i,m,t),dav(i,m,t)*rho(i,m,t));

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

Variables
zsp
y(i,k,t)
z(i,k,m,t)
a(i,k,t)
q(i,k,t)
b(i,k,t)
lvar(i,k,t);
Positive Variables y,z,a,q,b,lvar;

Equations
objsp
cons1_1sp(i,k,t)
cons1_2sp(i,k,t)
cons1_3sp(i,k,t)
cons2_1sp(i,k,t)
cons2_2sp(i,k,t)
cons3sp(i,k,t)
cons4_1sp(i,k,t)
cons4_2sp(i,k,t)
cons4_3sp(i,k,t)
cons5sp(i,k,t)
cons6sp(i,m,t);

objsp .. zsp =e= sum((i,k,t),c2(i)*y(i,k,t)+sum(m,(c3(i,k,m)-v1(i,t))*z(i,k,m,t)) + c4(i,t)*a(i,k,t) - v2(i,t)*q(i,k,t) + c5(i,k,t)*b(i,k,t)$(ord(t)<card(t)) + c6(i,k,t)*lvar(i,k,t)$(ord(t)<card(t)));
cons1_1sp(i,k,t)$(ord(t)=1) .. sum(j,x_star(i,j,k)) + a(i,k,t) =e= alpha(i)*y(i,k,t) + q(i,k,t) + b(i,k,t);
cons1_2sp(i,k,t)$(ord(t)>1 and ord(t)<card(t)) .. b(i,k,t-1) + a(i,k,t) =e= alpha(i)*y(i,k,t) + q(i,k,t) + b(i,k,t);
cons1_3sp(i,k,t)$(ord(t)=card(t)) .. b(i,k,t-1) + a(i,k,t) =e= alpha(i)*y(i,k,t) + q(i,k,t);
cons2_1sp(i,k,t)$(ord(t)=1) .. sum(j,x_star(i,j,k)) + a(i,k,t) =l= ca2(i,k);
cons2_2sp(i,k,t)$(ord(t)>1) .. b(i,k,t-1) + a(i,k,t) =l= ca2(i,k);
cons3sp(i,k,t)$(ord(t)<card(t)) .. lvar(i,k,t) =l= ca3(i,k);
cons4_1sp(i,k,t)$(ord(t)=1) .. y(i,k,t) =e= lvar(i,k,t) + sum(m,z(i,k,m,t));
cons4_2sp(i,k,t)$(ord(t)>1 and ord(t)<card(t)) .. lvar(i,k,t-1) + y(i,k,t) =e= lvar(i,k,t) + sum(m,z(i,k,m,t));
cons4_3sp(i,k,t)$(ord(t)=card(t)) .. lvar(i,k,t-1) + y(i,k,t) =e= sum(m,z(i,k,m,t));
cons5sp(i,k,t) .. y(i,k,t) =l= ca4av(i,k);
cons6sp(i,m,t) .. sum(k,z(i,k,m,t)) =l= dav(i,m,t);

Model subproblem
/objsp
cons1_1sp
cons1_2sp
cons1_3sp
cons2_1sp
cons2_2sp
cons3sp
cons4_1sp
cons4_2sp
cons4_3sp
cons5sp
cons6sp/;
option optcr=0;
option mip=cplex;
solve dualsubproblem maximizing zdsp using mip;
solve subproblem minimizing zsp using mip;
display zdsp.l, zsp.l, dualsubproblem.modelstat, subproblem.modelstat;

