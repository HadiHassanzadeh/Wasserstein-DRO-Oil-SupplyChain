from gurobipy import *
import numpy as np
import random

np.random.seed(42)
random.seed(42)

def power_function(x, power=0.3):
    return x ** power

I = range(3)
J = range(14)
K = range(6)
M = range(30)
T = range(12)
S = range(30)
N = range(30)
ITER = range(100)

x_max = 1
x_vals = np.linspace(0, x_max, len(N))
y_vals = power_function(x_vals)

f1_coeffs = {(i, j): random.randint(120, 250) for i in I for j in J}
ca1 = {(i, j): (random.randint(24, 30))/10 for i in I for j in J}
f1 = {(i, j, n): f1_coeffs[(i, j)] * float(y_vals[n-1])*ca1[(i,j)] for i in I for j in J for n in N}
g = {(i, j, n): ca1[(i, j)] * float(x_vals[n-1]) for i in I for j in J for n in N}

c1_price = {(i): ((random.randint(500, 520))/1000000) for i in I}
c1_distance = {(j,k): random.randint(300, 1000) for j in J for k in K}
c1 = {(i, j, k): c1_price[i]*c1_distance[j,k]*0.001 for i in I for j in J for k in K}
c2 = {(i): random.randint(69000, 71000)*0.000001 for i in I}
c3_price = {(i): (random.randint(800, 830))/1000000 for i in I for k in K for m in M}
c3_distance = {(k,m): random.randint(300, 1000) for k in K for m in M}
c3 = {(i,k,m): c3_price[i]*c3_distance[k,m]*0.001 for i in I for k in K for m in M}
c4 = {(i, t): random.randint(270000, 280000)*0.000001 for i in I for t in T}
c5 = {(i, k, t): random.randint(13000, 15000)*0.000001 for i in I for k in K for t in T}
c6 = {(i, k, t): random.randint(13000, 15000)*0.000001 for i in I for k in K for t in T}
v1 = {(i, t): random.randint(1000000, 1500000)*0.000001 for i in I for t in T}
v2 = {(i, t): random.randint(100000, 120000)*0.000001 for i in I for t in T}

ca2 = {(i, k): random.randint(95, 105) for i in I for k in K}
ca3 = {(i, k): random.randint(34, 38) for i in I for k in K}
base_ca4 = {(i, k): random.randint(32, 37) for i in I for k in K}
ca4 = {(i, k, s): max(1, int(np.random.normal(loc=base_ca4[i, k], scale=15))) for i in I for k in K for s in S}
ca4_upper = {(i, k): max(ca4[i,k,s] for s in S) for i in I for k in K}
ca4_lower = {(i, k): min(ca4[i,k,s] for s in S) for i in I for k in K}

base_dem = {(i, m, t): random.randint(1, 3) for i in I for m in M for t in T}
dem  = {(i, m, t, s): max(0, int(np.random.normal(loc=base_dem[i, m, t], scale=4))) for i in I for m in M for t in T for s in S}
dem_upper = {(i, m, t): max(dem[i,m,t,s] for s in S) for i in I for m in M for t in T}
dem_lower = {(i, m, t): min(dem[i,m,t,s] for s in S) for i in I for m in M for t in T}
alpha = {(i): 1/(random.randint(30, 35)/100) for i in I}

u = {(j, jp): 1 for j in J for jp in J}
p = {s: 1 / len(S) for s in S}
BigM = 10000

DM = Model("Deterministic Model")
DM.setParam('OutputFlag', 0)

delta = DM.addVars(J, vtype=GRB.BINARY, name="delta")
zeta = DM.addVars(I, J, N, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="zeta")
x = DM.addVars(I, J, K, lb=0, vtype=GRB.CONTINUOUS, name="x")
y = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="y")
z = DM.addVars(I, K, M, T, lb=0, vtype=GRB.CONTINUOUS, name="z")
a = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="a")
b = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="b")
l = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="l")
q = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="q")

obj_fn = (
    quicksum(f1[i,j,n]*zeta[i,j,n] for i in I for j in J for n in N)
    + quicksum(c1[i,j,k]*x[i,j,k] for i in I for j in J for k in K)
    + quicksum(c2[i]*y[i,k,t] + quicksum((c3[i,k,m]-v1[i,t])*z[i,k,m,t] for m in M) + c4[i,t]*a[i,k,t] - v2[i,t]*q[i,k,t] for i in I for k in K for t in T)
    + quicksum(c5[i,k,t]*b[i,k,t] + c6[i,k,t]*l[i,k,t] for i in I for k in K for t in T if t != len(T) - 1)
)

DM.setObjective(obj_fn, GRB.MINIMIZE)

for i in I:
    for j in J:
        DM.addSOS(GRB.SOS_TYPE2, [zeta[i, j, n] for n in N])


for i in I:
    for j in J:
        # --- Constraint 1 ---
        DM.addConstr(quicksum(x[i,j,k] for k in K) <= ca1[i,j])

        # --- Constraint 3 ---
        DM.addConstr(quicksum(x[i,j,k] for k in K) == quicksum(g[i,j,n]*zeta[i,j,n] for n in N))

        # --- Constraint 5 ---
        for k in K:
            DM.addConstr(x[i,j,k] <= BigM*delta[j])
        
        # --- Constraint 6 ---
        DM.addConstr(quicksum(zeta[i,j,n] for n in N) == 1)

# --- Constraint 2 ---
for i in I:
    for k in K:
        DM.addConstr(quicksum(x[i,j,k] for j in J) <= ca2[i,k])

# --- Constraint 4 ---
for j in J:
    for jp in J:
        DM.addConstr(delta[j] + delta[jp] <= 1 + u[j,jp])

for i in I:
    for k in K:
        for t in T:
            # --- Constraint 9 ---
            DM.addConstr(l[i,k,t] <= ca3[i,k])

            # --- Constraint 11 ---
            DM.addConstr(y[i,k,t] <= ca4[i,k])

            if t == 0:
                # --- Constraint 7.1 ---
                DM.addConstr(quicksum(x[i,j,k] for j in J) + a[i,k,t] == y[i,k,t]*alpha[i] + q[i,k,t] + b[i,k,t])

                # --- Constraint 8.1 ---
                DM.addConstr(quicksum(x[i,j,k] for j in J) + a[i,k,t] <= ca2[i,k])

                # --- Constraint 10.1 ---
                DM.addConstr(y[i,k,t] == l[i,k,t] + quicksum(z[i,k,m,t] for m in M))

            elif t == len(T) - 1:
                # --- Constraint 7.3 ---
                DM.addConstr(b[i,k,t-1] + a[i,k,t] == y[i,k,t]*alpha[i] + q[i,k,t])

                # --- Constraint 8.2 ---
                DM.addConstr(b[i,k,t-1] + a[i,k,t] <= ca2[i,k])

                # --- Constraint 10.3 ---
                DM.addConstr(l[i,k,t-1] + y[i,k,t] == quicksum(z[i,k,m,t] for m in M))

            else:
                # --- Constraint 7.2 ---
                DM.addConstr(b[i,k,t-1] + a[i,k,t] == y[i,k,t]*alpha[i] + q[i,k,t] + b[i,k,t])

                # --- Constraint 8.2 ---
                DM.addConstr(b[i,k,t-1] + a[i,k,t] <= ca2[i,k])

                # --- Constraint 10.2 ---
                DM.addConstr(l[i,k,t-1] + y[i,k,t] == l[i,k,t] + quicksum(z[i,k,m,t] for m in M))

# --- Constraint 12 ---
for i in I:
    for m in M:
        for t in T:
            DM.addConstr(quicksum(z[i,k,m,t] for k in K) <= dem[i,m,t])
                
DM.optimize()

# Check if the model has an optimal solution
if DM.status == GRB.OPTIMAL:
    print("🔹 Objective value:", DM.ObjVal)
else:
    print("❗ The model did not solve to optimality. Status code:", DM.status)

# print("\n🔹 Variable values (non-zero only):")
# for v in DM.getVars():
#     if v.X > 1e-6:  # Skip near-zero values
#         print(f"{v.VarName} = {v.X}")