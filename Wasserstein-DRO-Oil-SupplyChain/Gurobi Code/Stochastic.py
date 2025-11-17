from gurobipy import *
import numpy as np
import random

np.random.seed(42)
random.seed(42)

def g_function(x):
    Dalet = 0.5
    FinalValue = 1.0
    return FinalValue * (1 - np.exp(-Dalet * x))

I = range(2)     
J = range(5)    
K = range(3)  
M = range(15)  
T = range(12)    
S = range(20)   
N = range(60)


x_max = 60
x_vals = np.linspace(0, x_max, 61)
y_vals = g_function(x_vals)


f1_coeffs = {(i, j): random.randint(5, 10) for i in I for j in J}
g_coeffs = {(i, j): random.randint(10, 15) for i in I for j in J}
f1 = {(i, j, n): f1_coeffs[(i, j)] * float(y_vals[n - 1]) for i in I for j in J for n in N}
g = {(i, j, n): g_coeffs[(i, j)] * float(x_vals[n - 1]) for i in I for j in J for n in N}

c1 = {(i, j, k): random.randint(1, 5) for i in I for j in J for k in K}
c2 = {(i): random.randint(3, 7) for i in I}
c3 = {(i, k, m): random.randint(2, 6) for i in I for k in K for m in M}
c4 = {(i, t): random.randint(10, 15) for i in I for t in T}
c5 = {(i, k, t): random.randint(1, 3) for i in I for k in K for t in T}
c6 = {(i, k, t): random.randint(2, 4) for i in I for k in K for t in T}
v1 = {(i, t): random.randint(20, 30) for i in I for t in T}
v2 = {(i, t): random.randint(2, 8) for i in I for t in T}

ca1 = {(i, j): random.randint(40, 50) for i in I for j in J}
ca2 = {(i, k): random.randint(125, 150) for i in I for k in K}
ca3 = {(i, k): random.randint(175, 200) for i in I for k in K}
base_ca4 = {(i, k): random.randint(50, 75) for i in I for k in K}
ca4 = {(i, k, s): max(1, int(np.random.normal(loc=base_ca4[i, k], scale=15))) for i in I for k in K for s in S}

base_dem = {(i, m, t): random.randint(5, 10) for i in I for m in M for t in T}
dem  = {(i, m, t, s): max(0, int(np.random.normal(loc=base_dem[i, m, t], scale=4))) for i in I for m in M for t in T for s in S}
alpha = {(i): random.randint(2, 4) for i in I}

u = {(j, jp): 1 for j in J for jp in J}
pro = {s: 1 / len(S) for s in S}
BigM = 1000

SM = Model("Deterministic Model")
SM.setParam('OutputFlag', 0)

delta = SM.addVars(J, vtype=GRB.BINARY, name="delta")
zeta = SM.addVars(I, J, N, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="zeta")
x = SM.addVars(I, J, K, lb=0, vtype=GRB.CONTINUOUS, name="x")
y = SM.addVars(I, K, T, S, lb=0, vtype=GRB.CONTINUOUS, name="y")
z = SM.addVars(I, K, M, T, S, lb=0, vtype=GRB.CONTINUOUS, name="z")
a = SM.addVars(I, K, T, S, lb=0, vtype=GRB.CONTINUOUS, name="a")
b = SM.addVars(I, K, T, S, lb=0, vtype=GRB.CONTINUOUS, name="b")
l = SM.addVars(I, K, T, S, lb=0, vtype=GRB.CONTINUOUS, name="l")
q = SM.addVars(I, K, T, S, lb=0, vtype=GRB.CONTINUOUS, name="q")

obj_fn = (
    quicksum(f1[i,j,n]*zeta[i,j,n] for i in I for j in J for n in N)
    + quicksum(c1[i,j,k]*x[i,j,k] for i in I for j in J for k in K)
    + quicksum(pro[s]*(
        quicksum(c2[i]*y[i,k,t,s] + quicksum((c3[i,k,m]-v1[i,t])*z[i,k,m,t,s] for m in M) + c4[i,t]*a[i,k,t,s] - v2[i,t]*q[i,k,t,s] for i in I for k in K for t in T)
        + quicksum(c5[i,k,t]*b[i,k,t,s] + c6[i,k,t]*l[i,k,t,s] for i in I for k in K for t in T if t != len(T) - 1)
    ) for s in S)
)

SM.setObjective(obj_fn, GRB.MINIMIZE)

for i in I:
    for j in J:
        SM.addSOS(GRB.SOS_TYPE2, [zeta[i, j, n] for n in N])


for i in I:
    for j in J:
        # --- Constraint 1 ---
        SM.addConstr(quicksum(x[i,j,k] for k in K) <= ca1[i,j])

        # --- Constraint 3 ---
        SM.addConstr(quicksum(x[i,j,k] for k in K) == quicksum(g[i,j,n]*zeta[i,j,n] for n in N))

        # --- Constraint 5 ---
        for k in K:
            SM.addConstr(x[i,j,k] <= BigM*delta[j])
        
        # --- Constraint 6 ---
        SM.addConstr(quicksum(zeta[i,j,n] for n in N) == 1)

# --- Constraint 2 ---
for i in I:
    for k in K:
        SM.addConstr(quicksum(x[i,j,k] for j in J) <= ca2[i,k])

# --- Constraint 4 ---
for j in J:
    for jp in J:
        SM.addConstr(delta[j] + delta[jp] <= 1 + u[j,jp])

for i in I:
    for k in K:
        for t in T:
            for s in S:
                # --- Constraint 9 ---
                SM.addConstr(l[i,k,t,s] <= ca3[i,k])

                # --- Constraint 11 ---
                SM.addConstr(y[i,k,t,s] <= ca4[i,k,s])

                if t == 0:
                    # --- Constraint 7.1 ---
                    SM.addConstr(quicksum(x[i,j,k] for j in J) + a[i,k,t,s] == y[i,k,t,s]*alpha[i] + q[i,k,t,s] + b[i,k,t,s])

                    # --- Constraint 8.1 ---
                    SM.addConstr(quicksum(x[i,j,k] for j in J) + a[i,k,t,s] <= ca2[i,k])

                    # --- Constraint 10.1 ---
                    SM.addConstr(y[i,k,t,s] == l[i,k,t,s] + quicksum(z[i,k,m,t,s] for m in M))

                elif t == len(T) - 1:
                    # --- Constraint 7.3 ---
                    SM.addConstr(b[i,k,t-1,s] + a[i,k,t,s] == y[i,k,t,s]*alpha[i] + q[i,k,t,s])

                    # --- Constraint 8.2 ---
                    SM.addConstr(b[i,k,t-1,s] + a[i,k,t,s] <= ca2[i,k])

                    # --- Constraint 10.3 ---
                    SM.addConstr(l[i,k,t-1,s] + y[i,k,t,s] == quicksum(z[i,k,m,t,s] for m in M))

                else:
                    # --- Constraint 7.2 ---
                    SM.addConstr(b[i,k,t-1,s] + a[i,k,t,s] == y[i,k,t,s]*alpha[i] + q[i,k,t,s] + b[i,k,t,s])

                    # --- Constraint 8.2 ---
                    SM.addConstr(b[i,k,t-1,s] + a[i,k,t,s] <= ca2[i,k])

                    # --- Constraint 10.2 ---
                    SM.addConstr(l[i,k,t-1,s] + y[i,k,t,s] == l[i,k,t,s] + quicksum(z[i,k,m,t,s] for m in M))

# --- Constraint 12 ---
for i in I:
    for m in M:
        for t in T:
            for s in S:
                SM.addConstr(quicksum(z[i,k,m,t,s] for k in K) <= dem[i,m,t,s])
                
SM.optimize()

# Check if the model has an optimal solution
if SM.status == GRB.OPTIMAL:
    print("🔹 Objective value:", SM.ObjVal)
else:
    print("❗ The model did not solve to optimality. Status code:", SM.status)

print("\n🔹 Variable values (non-zero only):")
for v in SM.getVars():
    if v.X > 1e-6:  # Skip near-zero values
        print(f"{v.VarName} = {v.X}")
