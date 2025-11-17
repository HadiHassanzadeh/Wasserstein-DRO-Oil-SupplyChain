from gurobipy import *
import numpy as np
import random

def g_function(x):
    Dalet = 0.5
    FinalValue = 1.0
    return FinalValue * (1 - np.exp(-Dalet * x))

I = range(2)     
J = range(5)    
K = range(3)  
M = range(15)  
T = range(12)    
S = range(100)   
N = range(60) 


x_max = 60
x_vals = np.linspace(0, x_max, 61)
y_vals = g_function(x_vals)


f1_coeffs = {(i, j, s): random.randint(5, 10) for i in I for j in J for s in S}
g_coeffs = {(i, j, s): random.randint(10, 15) for i in I for j in J for s in S}
f1 = {(i, j, n, s): f1_coeffs[(i, j, s)] * float(y_vals[n - 1]) for i in I for j in J for n in N for s in S}
g = {(i, j, n, s): g_coeffs[(i, j, s)] * float(x_vals[n - 1]) for i in I for j in J for n in N for s in S}

c1 = {(i, j, k, s): random.randint(1, 5) for i in I for j in J for k in K for s in S}
c2 = {(i, s): random.randint(3, 7) for i in I for s in S}
c3 = {(i, k, m, s): random.randint(2, 6) for i in I for k in K for m in M for s in S}
c4 = {(i, t, s): random.randint(10, 15) for i in I for t in T for s in S}
c5 = {(i, k, t, s): random.randint(1, 3) for i in I for k in K for t in T for s in S}
c6 = {(i, k, t, s): random.randint(2, 4) for i in I for k in K for t in T for s in S}
v1 = {(i, t, s): random.randint(20, 30) for i in I for t in T for s in S}
v2 = {(i, t, s): random.randint(2, 8) for i in I for t in T for s in S}

ca1 = {(i, j, s): random.randint(40, 50) for i in I for j in J for s in S}
ca2 = {(i, k, s): random.randint(125, 150) for i in I for k in K for s in S}
ca3 = {(i, k, s): random.randint(175, 200) for i in I for k in K for s in S}
ca4 = {(i, k, s): random.randint(50, 75) for i in I for k in K for s in S}

dem = {(i, m, t, s): random.randint(3, 10) for i in I for m in M for t in T for s in S}
alpha = {(i, s): random.randint(2, 4) for i in I for s in S}

u = {(j, jp, s): 1 for j in J for jp in J for s in S}
p = {s: 1 / len(S) for s in S}
BigM = 1000
x_star = {(i, j, k): 1 for i in I for j in J for k in K}

def Primal(sa):
    PM = Model("Primal Model")
    PM.setParam('OutputFlag', 0)

    y = PM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="y")
    z = PM.addVars(I, K, M, T, lb=0, vtype=GRB.CONTINUOUS, name="z")
    a = PM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="a")
    b = PM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="b")
    l = PM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="l")
    q = PM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="q")

    obj_fn = quicksum(
        c2[i, sa] * y[i, k, t] +
        quicksum((c3[i, k, m, sa] - v1[i, t, sa]) * z[i, k, m, t] for m in M) +
        c4[i, t, sa] * a[i, k, t] -
        v2[i, t, sa] * q[i, k, t] +
        c5[i, k, t, sa] * b[i, k, t] +
        c6[i, k, t, sa] * l[i, k, t]
        for i in I for k in K for t in T
    )

    PM.setObjective(obj_fn, GRB.MINIMIZE)

    for i in I:
        for k in K:
            for t in T:
                # --- Constraint 1 ---
                if t == 0:
                    lhs_c1 = quicksum(x_star[i, j, k] for j in J) + a[i, k, t]
                    rhs_c1 = y[i, k, t] * alpha[i, sa] + q[i, k, t] + b[i, k, t]
                elif t == len(T) - 1:
                    lhs_c1 = b[i, k, t - 1] + a[i, k, t]
                    rhs_c1 = y[i, k, t] * alpha[i, sa] + q[i, k, t]
                else:
                    lhs_c1 = b[i, k, t - 1] + a[i, k, t]
                    rhs_c1 = y[i, k, t] * alpha[i, sa] + q[i, k, t] + b[i, k, t]
                PM.addConstr(lhs_c1 == rhs_c1, name=f"c1_{i}_{k}_{t}")

                # --- Constraint 2 ---
                if t == 0:
                    lhs_c2 = quicksum(x_star[i, j, k] for j in J) + a[i, k, t]
                else:
                    lhs_c2 = b[i, k, t - 1] + a[i, k, t]
                PM.addConstr(lhs_c2 <= ca2[i, k, sa], name=f"c2_{i}_{k}_{t}")

                # --- Constraint 3 ---
                PM.addConstr(l[i, k, t] <= ca3[i, k, sa])

                # --- Constraint 4 ---
                if t == 0:
                    lhs_c4 = y[i, k, t]
                    rhs_c4 = l[i, k, t] + quicksum(z[i, k, m, t] for m in M)
                elif t == len(T) - 1:
                    lhs_c4 = l[i, k, t - 1] + y[i, k, t]
                    rhs_c4 = quicksum(z[i, k, m, t] for m in M)
                else:
                    lhs_c4 = l[i, k, t - 1] + y[i, k, t]
                    rhs_c4 = l[i, k, t] + quicksum(z[i, k, m, t] for m in M)
                PM.addConstr(lhs_c4 == rhs_c4, name=f"c4_{i}_{k}_{t}")

                # --- Constraint 5 ---
                PM.addConstr(y[i, k, t] <= ca4[i, k, sa])

    # --- Constraint 6 ---
    for i in I:
        for m in M:
            for t in T:
                PM.addConstr(quicksum(z[i, k, m, t] for k in K) <= dem[i, m, t, sa], name=f"c6_{i}_{m}_{t}")
    
    PM.optimize()
    
    Primal_Obj_Value = None
    if PM.status == GRB.OPTIMAL:
        Primal_Obj_Value = PM.ObjVal
    
    return Primal_Obj_Value

def Dual(sa):
    DM = Model("Dual Model")
    DM.setParam('OutputFlag', 0)

    gamma = DM.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="gamma")
    mu = DM.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="mu")
    eta = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="eta")
    nu = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="nu")
    rho = DM.addVars(I, M, T, lb=0, vtype=GRB.CONTINUOUS, name="rho")
    pi = DM.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="pi")

    obj_fn = (-1*quicksum(ca4[i,k,sa]*nu[i,k,t] for i in I for k in K for t in T) 
            - quicksum((ca2[i,k,sa]-sum(x_star[i,j,k] for j in J))*eta[i,k,0] for i in I for k in K) 
            - quicksum(ca2[i,k,sa]*eta[i,k,t] for i in I for k in K for t in T if t != 0) 
            - quicksum(x_star[i,j,k]*gamma[i,k,0] for i in I for j in J for k in K) 
            - quicksum(ca3[i,k,sa]*pi[i,k,t] for i in I for k in K for t in T if t != len(T) - 1) 
            - quicksum(dem[i,m,t,sa]*rho[i,m,t] for i in I for m in M for t in T))
    
    DM.setObjective(obj_fn, GRB.MAXIMIZE)

    for i in I:
        for k in K:
            for t in T:
                # --- Constraint 1 ---
                DM.addConstr(-1*alpha[i,sa]*gamma[i,k,t] + mu[i,k,t] - nu[i,k,t] <= c2[i,sa], name=f'c1_{i}_{k}_{t}')

                # --- Constraint 2 ---
                for m in M:
                    DM.addConstr(-1*mu[i,k,t]- rho[i,m,t] <= c3[i,k,m,sa] - v1[i,t,sa], name=f'c2_{i}_{k}_{m}_{t}')

                # --- Constraint 3 ---
                DM.addConstr(gamma[i,k,t] + eta[i,k,t] <= c4[i,t,sa], name=f'c3_{i}_{k}_{t}')

                # --- Constraint 4 ---
                DM.addConstr(gamma[i,k,t] >= v2[i,t,sa], name=f'c4_{i}_{k}_{t}')

                if t != len(T) - 1:
                    # --- Constraint 5 ---
                    DM.addConstr(-1*gamma[i,k,t] + gamma[i,k,t+1] + eta[i,k,t+1] <= c5[i,k,t,sa], name=f'c5_{i}_{k}_{t}')

                    # --- Constraint 6 ---
                    DM.addConstr(-1*mu[i,k,t] + mu[i,k,t+1] + pi[i,k,t] <= c6[i,k,t,sa], name=f'c6_{i}_{k}_{t}')
    
    DM.optimize()
    
    Dual_Obj_Value = None
    if DM.status == GRB.OPTIMAL:
        Dual_Obj_Value = DM.ObjVal
    
    return Dual_Obj_Value

for sa in S:
    print(10*"=")
    p = Primal(sa)
    d = Dual(sa)
    print(f"Sampel Number: {sa + 1}")
    print(f"Primal Obj: {p}")
    print(f"Dual Obj: {d}")
    if p == d:
        print("✅ OK")
    else:
        print("❌ NOT OK") 