sets
i /1*3/
j /1*5/
k /1*3/
m /1*15/
t /1*12/
sa /1*10/
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

*////////////////////////// Stochastic Model ///////////////////////
Variables
zs
y(i,k,t,sa)
z(i,k,m,t,sa)
a(i,k,t,sa)
q(i,k,t,sa)
b(i,k,t,sa)
lvar(i,k,t,sa)
zeta(i,j,n)
x(i,j,k)
delta(j);
Free Variables zs;
Sos2 Variables zeta;
Positive Variables x, y, z, a, q, b, lvar;
Binary Variables delta;
Equations
obj
cons1(i,j)
cons2(i,k)
cons3(i,j)
cons4(j,jp)
cons5(i,j,k)
cons6(i,j)
cons7_1(i,k,t,sa)
cons7_2(i,k,t,sa)
cons7_3(i,k,t,sa)
cons8_1(i,k,t,sa)
cons8_2(i,k,t,sa)
cons9(i,k,t,sa)
cons10_1(i,k,t,sa)
cons10_2(i,k,t,sa)
cons10_3(i,k,t,sa)
cons11(i,k,t,sa)
cons12(i,m,t,sa);

obj .. zs =e= sum((i,j,n),f1(i,j,n)*zeta(i,j,n)) + sum((i,j,k),c1(i,j,k)*x(i,j,k)) + sum(sa,sampel_size*(sum((i,k,t),c2(i)*y(i,k,t,sa)) + sum((i,k,t,m),(c3(i,k,m)-v1(i,t))*z(i,k,m,t,sa)) + sum((i,k,t),c4(i,t)*a(i,k,t,sa)) - sum((i,k,t),v2(i,t)*q(i,k,t,sa)) + sum((i,k,t)$(ord(t)<card(t)),c5(i,k,t)*b(i,k,t,sa)) + sum((i,k,t)$(ord(t)<card(t)),c6(i,k,t)*lvar(i,k,t,sa))));

cons1(i,j) .. sum(k,x(i,j,k)) =l= ca1(i,j);
cons2(i,k) .. sum(j,x(i,j,k)) =l= ca2(i,k);
cons3(i,j) .. sum(k,x(i,j,k)) =e= sum(n,g(i,j,n)*zeta(i,j,n));
cons4(j,jp) .. delta(j) + delta(jp) =l= 1 + u(j,jp);
cons5(i,j,k) .. x(i,j,k) =l= BigM*delta(j);
cons6(i,j) .. sum(n,zeta(i,j,n)) =e= 1;
cons7_1(i,k,t,sa)$(ord(t)=1) .. sum(j,x(i,j,k)) + a(i,k,t,sa) =e= alpha(i)*y(i,k,t,sa) + q(i,k,t,sa) + b(i,k,t,sa);
cons7_2(i,k,t,sa)$(ord(t)>1 and ord(t)<card(t)) .. b(i,k,t-1,sa) + a(i,k,t,sa) =e= alpha(i)*y(i,k,t,sa) + q(i,k,t,sa) + b(i,k,t,sa);
cons7_3(i,k,t,sa)$(ord(t)=card(t)) .. b(i,k,t-1,sa) + a(i,k,t,sa) =e= alpha(i)*y(i,k,t,sa) + q(i,k,t,sa);
cons8_1(i,k,t,sa)$(ord(t)=1) .. sum(j,x(i,j,k)) + a(i,k,t,sa) =l= ca2(i,k);
cons8_2(i,k,t,sa)$(ord(t)>1) .. b(i,k,t-1,sa) + a(i,k,t,sa) =l= ca2(i,k);
cons9(i,k,t,sa)$(ord(t)<card(t)) .. lvar(i,k,t,sa) =l= ca3(i,k);
cons10_1(i,k,t,sa)$(ord(t)=1) .. y(i,k,t,sa) =e= lvar(i,k,t,sa) + sum(m,z(i,k,m,t,sa));
cons10_2(i,k,t,sa)$(ord(t)>1 and ord(t)<card(t)) .. lvar(i,k,t-1,sa) + y(i,k,t,sa) =e= lvar(i,k,t,sa) + sum(m,z(i,k,m,t,sa));
cons10_3(i,k,t,sa)$(ord(t)=card(t)) .. lvar(i,k,t-1,sa) + y(i,k,t,sa) =e= sum(m,z(i,k,m,t,sa));
cons11(i,k,t,sa) .. y(i,k,t,sa) =l= ca4(i,k,sa);
cons12(i,m,t,sa) .. sum(k,z(i,k,m,t,sa)) =l= d(i,m,t,sa);

Model stochastic
/obj,
cons1,
cons2,
cons3,
cons4,
cons5,
cons6,
cons7_1,
cons7_2,
cons7_3,
cons8_1,
cons8_2,
cons9,
cons10_1,
cons10_2,
cons10_3,
cons11,
cons12/;

Scalars startTime, endTime, executionTime;
startTime = jnow;
option optcr=0;
option mip=cplex;
solve stochastic minimizing zs using mip;
endTime = jnow;
executionTime = (startTime-endTime)*86400;
display zs.l, x.l, zeta.l, delta.l, executionTime;
