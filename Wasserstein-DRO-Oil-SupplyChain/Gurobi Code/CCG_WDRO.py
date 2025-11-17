from gurobipy import *
import numpy as np
import random
import time

np.random.seed(42)
random.seed(42)

def power_function(x, power=0.3):
    return x ** power

def g_function(x):
    Dalet = 0.5
    FinalValue = 1.0
    return FinalValue * (1 - np.exp(-Dalet * x))

theta = 700
I = range(3)
J = range(15)
K = range(5)
M = range(30)
S = range(10)
SST = range(500)
T = range(12)
N = range(60)
ITER = range(500)
# S_Total = set(range(100))
# SST = set(range(0, 20))
# S = S_Total - SST


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
#ca4_pool = {(i, k, s): max(1, int(np.random.normal(loc=base_ca4[i, k], scale=15))) for i in I for k in K for s in S_Total}
#ca4 = {(i, k, s): max(1, int(np.random.exponential(scale=base_ca4[i, k]))) for i in I for k in K for s in S}


base_dem = {(i, m, t): 7 for i in I for m in M for t in T}
#dem_pool  = {(i, m, t, s): max(0, int(np.random.normal(loc=base_dem[i, m, t], scale=4))) for i in I for m in M for t in T for s in S_Total}
#dem  = {(i, m, t, s): max(0, int(np.random.exponential(scale=base_dem[i, m, t]))) for i in I for m in M for t in T for s in S}


alpha = {(i): random.randint(2, 4) for i in I}

ca4 = {
    (i, k, s): max(1, int(np.random.normal(loc=base_ca4[i, k], scale=20)))
    for i in I for k in K for s in S
}

dem = {
    (i, m, t, s): max(0, int(np.random.normal(loc=base_dem[i, m, t], scale=20)))
    for i in I for m in M for t in T for s in S
}

ca4_test = {
    (i, k, s): max (1, np.random.poisson(lam=base_ca4[i,k]))
    for i in I for k in K for s in SST}

dem_test = {
    (i, m, t, s): max(0, np.random.poisson(lam=base_dem[i, m, t]))
    for i in I for m in M for t in T for s in SST
}

ca4_upper = {(i, k): max(ca4[i,k,s] for s in S) for i in I for k in K}
ca4_lower = {(i, k): min(ca4[i,k,s] for s in S) for i in I for k in K}
dem_upper = {(i, m, t): max(dem[i,m,t,s] for s in S) for i in I for m in M for t in T}
dem_lower = {(i, m, t): min(dem[i,m,t,s] for s in S) for i in I for m in M for t in T}

u = {(j, jp): 1 for j in J for jp in J}
pro = {s: 1 / len(S) for s in S}
BigM = 1000

def Dual_Sub_Problem(sample,
                     I, J, K, T, M,
                     ca2, ca3, ca4, c2, c3, c4, c5, c6,
                     v1, v2, alpha, dem, x_star, landa_star, ca4_upper, ca4_lower, dem_upper, dem_lower):
    DSP = Model("Dual Sub Problem")
    DSP.setParam('OutputFlag', 0)

    gamma = DSP.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="gamma")
    mu = DSP.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="mu")
    eta = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="eta")
    nu = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="nu")
    rho = DSP.addVars(I, M, T, lb=0, vtype=GRB.CONTINUOUS, name="rho")
    pi = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="pi")
    omega_plus = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name='omega_plus')
    omega_minus = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name='omega_minus')
    OMEGA_plus = DSP.addVars(I, M, T, lb=0, vtype=GRB.CONTINUOUS, name="OMEGA_plus")
    OMEGA_minus = DSP.addVars(I, M, T, lb=0, vtype=GRB.CONTINUOUS, name="OMEGA_minus")
    chi_plus = DSP.addVars(I, K, T, vtype=GRB.BINARY, name="chi_plus")
    chi_minus = DSP.addVars(I, K, T, vtype=GRB.BINARY, name="chi_minus")
    GAMMA_plus = DSP.addVars(I, M, T, vtype=GRB.BINARY, name="GAMMA_plus")
    GAMMA_minus = DSP.addVars(I, M, T, vtype=GRB.BINARY, name="GAMMA_minus")



    obj_fn = (-1*quicksum((ca4_upper[i,k] - ca4[i,k,sample])*omega_plus[i,k,t] + (ca4_lower[i,k] - ca4[i,k,sample])*omega_minus[i,k,t] + ca4[i,k,sample]*nu[i,k,t] for i in I for k in K for t in T) 
            - quicksum((ca2[i,k]-sum(x_star[i,j,k] for j in J))*eta[i,k,0] for i in I for k in K) 
            - quicksum(ca2[i,k]*eta[i,k,t] for i in I for k in K for t in T if t != 0) 
            - quicksum(x_star[i,j,k]*gamma[i,k,0] for i in I for j in J for k in K) 
            - quicksum(ca3[i,k]*pi[i,k,t] for i in I for k in K for t in T if t != len(T) - 1) 
            - quicksum((dem_upper[i,m,t] - dem[i,m,t,sample])*OMEGA_plus[i,m,t] + (dem_lower[i,m,t] - dem[i,m,t,sample])*OMEGA_minus[i,m,t] + dem[i,m,t,sample]*rho[i,m,t] for i in I for m in M for t in T)
            - quicksum(landa_star*((ca4_upper[i,k] - ca4[i,k,sample])*chi_plus[i,k,t] + (ca4[i,k,sample] - ca4_lower[i,k])*chi_minus[i,k,t]) for i in I for k in K for t in T)
            - quicksum(landa_star*((dem_upper[i,m,t] - dem[i,m,t,sample])*GAMMA_plus[i,m,t] + (dem[i,m,t,sample] - dem_lower[i,m,t])*GAMMA_minus[i,m,t]) for i in I for m in M for t in T))
    
    DSP.setObjective(obj_fn, GRB.MAXIMIZE)

    for i in I:
        for k in K:
            for t in T:
                # --- Constraint 1 ---
                DSP.addConstr(-1*alpha[i]*gamma[i,k,t] + mu[i,k,t] - nu[i,k,t] <= c2[i])

                # --- Constraint 2 ---
                for m in M:
                    DSP.addConstr(-1*mu[i,k,t]- rho[i,m,t] <= c3[i,k,m] - v1[i,t])

                # --- Constraint 3 ---
                DSP.addConstr(gamma[i,k,t] + eta[i,k,t] <= c4[i,t])

                # --- Constraint 4 ---
                DSP.addConstr(gamma[i,k,t] >= v2[i,t])

                if t != len(T) - 1:
                    # --- Constraint 5 ---
                    DSP.addConstr(-1*gamma[i,k,t] + gamma[i,k,t+1] + eta[i,k,t+1] <= c5[i,k,t])

                    # --- Constraint 6 ---
                    DSP.addConstr(-1*mu[i,k,t] + mu[i,k,t+1] + pi[i,k,t] <= c6[i,k,t])
    
    for i in I:
        for k in K:
            for t in T:
                # --- Constraint 7 ---
                DSP.addConstr(omega_plus[i,k,t] <= BigM*chi_plus[i,k,t])

                # --- Constraint 8 ---
                DSP.addConstr(omega_plus[i,k,t] <= nu[i,k,t])

                # --- Constraint 9 ---
                DSP.addConstr(omega_plus[i,k,t] >= nu[i,k,t] - BigM*(1-chi_plus[i,k,t]))

                # --- Constraint 10 ---
                DSP.addConstr(omega_minus[i,k,t] <= BigM*chi_minus[i,k,t])

                # --- Constraint 11 ---
                DSP.addConstr(omega_minus[i,k,t] <= nu[i,k,t])

                # --- Constraint 12 ---
                DSP.addConstr(omega_minus[i,k,t] >= nu[i,k,t] - BigM*(1-chi_minus[i,k,t]))

                # --- Constraint 13 ---
                DSP.addConstr(chi_plus[i,k,t] + chi_minus[i,k,t] <= 1)
    
    for i in I:
        for m in M:
            for t in T:
                # --- Constraint 14 ---
                DSP.addConstr(OMEGA_plus[i,m,t] <= BigM*GAMMA_plus[i,m,t])

                # --- Constraint 15 ---
                DSP.addConstr(OMEGA_plus[i,m,t] <= rho[i,m,t])

                # --- Constraint 16 ---
                DSP.addConstr(OMEGA_plus[i,m,t] >= rho[i,m,t] - BigM*(1-GAMMA_plus[i,m,t]))

                # --- Constraint 17 ---
                DSP.addConstr(OMEGA_minus[i,m,t] <= BigM*GAMMA_minus[i,m,t])

                # --- Constraint 18 ---
                DSP.addConstr(OMEGA_minus[i,m,t] <= rho[i,m,t])

                # --- Constraint 19 ---
                DSP.addConstr(OMEGA_minus[i,m,t] >= rho[i,m,t] - BigM*(1-GAMMA_minus[i,m,t]))

                # --- Constraint 20 ---
                DSP.addConstr(GAMMA_plus[i,m,t] + GAMMA_minus[i,m,t] <= 1)

    
    DSP.optimize()


    if DSP.status == GRB.OPTIMAL:

        DSP_Obj_Value = DSP.ObjVal
        chi_plus_sol = {(i,k,t): chi_plus[i,k,t].X for i in I for k in K for t in T}
        chi_minus_sol = {(i,k,t): chi_minus[i,k,t].X for i in I for k in K for t in T}
        GAMMA_plus_sol = {(i,m,t): GAMMA_plus[i,m,t].X for i in I for m in M for t in T}
        GAMMA_minus_sol = {(i,m,t): GAMMA_minus[i,m,t].X for i in I for m in M for t in T}
    

        return DSP_Obj_Value, chi_plus_sol, chi_minus_sol, GAMMA_plus_sol, GAMMA_minus_sol
    else:
        print("❌ DSP ERROR ❌")  
        return None, None, None, None, None

def Restricted_Master_Problem(iter,
                              I, J, K, N, S, T, M,
                              f1, c1, g, ca1, ca2, ca3, ca4, 
                              u, BigM, theta, ca4_upper, ca4_lower, dem_upper, dem_lower,
                              chi_plus_star, chi_minus_star, GAMMA_plus_star, GAMMA_minus_star):

    IT = range(iter+1)
    RMP = Model("Restricted Master Problem")
    RMP.setParam('OutputFlag', 0)

    delta = RMP.addVars(J, vtype=GRB.BINARY, name="delta")
    zeta = RMP.addVars(I, J, N, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="zeta")
    x = RMP.addVars(I, J, K, lb=0, vtype=GRB.CONTINUOUS, name="x")
    landa = RMP.addVar(lb=0, vtype=GRB.CONTINUOUS, name="landa")
    y = RMP.addVars(I, K, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="y")
    z = RMP.addVars(I, K, M, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="z")
    a = RMP.addVars(I, K, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="a")
    b = RMP.addVars(I, K, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="b")
    l = RMP.addVars(I, K, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="l")
    q = RMP.addVars(I, K, T, ITER, S, lb=0, vtype=GRB.CONTINUOUS, name="q")
    phi_cut = RMP.addVars(S, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="phi_cut")

    obj_fn = (
    quicksum(f1[i,j,n]*zeta[i,j,n] for i in I for j in J for n in N)
    + quicksum(c1[i,j,k]*x[i,j,k] for i in I for j in J for k in K)
    + landa*theta
    + (1/len(S))*quicksum(phi_cut[s] for s in S))

    RMP.setObjective(obj_fn, GRB.MINIMIZE)

    for i in I:
        for j in J:
            RMP.addSOS(GRB.SOS_TYPE2, [zeta[i, j, n] for n in N])

    for i in I:
        for j in J:
            # --- Constraint 1 FS ---
            RMP.addConstr(quicksum(x[i,j,k] for k in K) <= ca1[i,j])

            # --- Constraint 3  FS ---
            RMP.addConstr(quicksum(x[i,j,k] for k in K) == quicksum(g[i,j,n]*zeta[i,j,n] for n in N))

            # --- Constraint 5 FS ---
            for k in K:
                RMP.addConstr(x[i,j,k] <= BigM*delta[j])
            
            # --- Constraint 6 FS ---
            RMP.addConstr(quicksum(zeta[i,j,n] for n in N) == 1)

    # --- Constraint 2 FS ---
    for i in I:
        for k in K:
            RMP.addConstr(quicksum(x[i,j,k] for j in J) <= ca2[i,k])

    # --- Constraint 4 FS ---
    for j in J:
        for jp in J:
            RMP.addConstr(delta[j] + delta[jp] <= 1 + u[j,jp])

    # --- Optimality Cut ---
    for iter in IT:
        for sample in S:
            lhs = phi_cut[sample]
            rhs = (quicksum(
                    c2[i] * y[i, k, t, iter, sample] +
                    quicksum((c3[i, k, m] - v1[i, t]) * z[i, k, m, t, iter, sample] for m in M) +
                    c4[i, t] * a[i, k, t, iter, sample] -
                    v2[i, t] * q[i, k, t, iter, sample] +
                    c5[i, k, t] * b[i, k, t, iter, sample] +
                    c6[i, k, t] * l[i, k, t, iter, sample]
                    for i in I for k in K for t in T)
                    - quicksum(landa*((ca4_upper[i,k] - ca4[i,k,sample])*chi_plus_star[i,k,t,iter,sample] + (ca4[i,k,sample] - ca4_lower[i,k])*chi_minus_star[i,k,t,iter,sample]) for i in I for k in K for t in T)
                    - quicksum(landa*((dem_upper[i,m,t] - dem[i,m,t,sample])*GAMMA_plus_star[i,m,t,iter,sample] + (dem[i,m,t,sample] - dem_lower[i,m,t])*GAMMA_minus_star[i,m,t,iter,sample]) for i in I for m in M for t in T))
            RMP.addConstr(lhs >= rhs)
    
    for iter in ITER:
        for sample in S:
            for i in I:
                for k in K:
                    for t in T:
                        # --- Constraint 1 SS ---
                        if t == 0:
                            lhs_c1 = quicksum(x[i, j, k] for j in J) + a[i, k, t, iter, sample]
                            rhs_c1 = y[i, k, t, iter, sample] * alpha[i] + q[i, k, t, iter, sample] + b[i, k, t, iter, sample]
                        elif t == len(T) - 1:
                            lhs_c1 = b[i, k, t - 1, iter, sample] + a[i, k, t, iter, sample]
                            rhs_c1 = y[i, k, t, iter, sample] * alpha[i] + q[i, k, t, iter, sample]
                        else:
                            lhs_c1 = b[i, k, t - 1, iter, sample] + a[i, k, t, iter, sample]
                            rhs_c1 = y[i, k, t, iter, sample] * alpha[i] + q[i, k, t, iter, sample] + b[i, k, t, iter, sample]
                        RMP.addConstr(lhs_c1 == rhs_c1)

                        # --- Constraint 2 SS ---
                        if t == 0:
                            lhs_c2 = quicksum(x[i, j, k] for j in J) + a[i, k, t, iter, sample]
                        else:
                            lhs_c2 = b[i, k, t - 1, iter, sample] + a[i, k, t, iter, sample]
                        RMP.addConstr(lhs_c2 <= ca2[i, k])

                        # --- Constraint 3 SS ---
                        RMP.addConstr(l[i, k, t, iter, sample] <= ca3[i, k])

                        # --- Constraint 4 SS ---
                        if t == 0:
                            lhs_c4 = y[i, k, t, iter, sample]
                            rhs_c4 = l[i, k, t, iter, sample] + quicksum(z[i, k, m, t, iter, sample] for m in M)
                        elif t == len(T) - 1:
                            lhs_c4 = l[i, k, t - 1, iter, sample] + y[i, k, t, iter, sample]
                            rhs_c4 = quicksum(z[i, k, m, t, iter, sample] for m in M)
                        else:
                            lhs_c4 = l[i, k, t - 1, iter, sample] + y[i, k, t, iter, sample]
                            rhs_c4 = l[i, k, t, iter, sample] + quicksum(z[i, k, m, t, iter, sample] for m in M)
                        RMP.addConstr(lhs_c4 == rhs_c4)

                        # --- Constraint 5 SS ---
                        RMP.addConstr(y[i, k, t, iter, sample] <= (ca4[i, k, sample] + (ca4_upper[i,k] - ca4[i,k,sample])*chi_plus_star[i,k,t,iter,sample] + (ca4[i,k,sample] - ca4_lower[i,k])*chi_minus_star[i,k,t,iter,sample]))  

            # --- Constraint 6 SS ---
            for i in I:
                for m in M:
                    for t in T:
                        RMP.addConstr(quicksum(z[i, k, m, t, iter, sample] for k in K) <= (dem[i, m, t, sample] + (dem_upper[i,m,t] - dem[i,m,t,sample])*GAMMA_plus_star[i,m,t,iter,sample] + (dem[i,m,t,sample] - dem_lower[i,m,t])*GAMMA_minus_star[i,m,t,iter,sample]))
            
    # --- Feasibility Cut ---
    RMP.addConstr(landa <= 10000)
    
    RMP.optimize()

    if RMP.status == GRB.OPTIMAL:
        RMP_Obj_Value = RMP.ObjVal

        x_sol = {(i,j,k): x[i,j,k].X for i in I for j in J for k in K}
        zeta_sol = {(i,j,n): zeta[i,j,n].X for i in I for j in J for n in N}
        landa_sol = landa.X
        
        return RMP_Obj_Value, x_sol, zeta_sol, landa_sol

    else:
        print("❌ RMP ERROR ❌")  
        return None, None, None, None, None

Previous_UB = 0


landa_star = 0
zeta_star = np.zeros((len(I), len(J), len(N)))
x_star = np.zeros((len(I), len(J), len(K)))

chi_plus_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
chi_minus_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
GAMMA_plus_star = np.zeros((len(I), len(M), len(T), len(ITER), len(S)))
GAMMA_minus_star = np.zeros((len(I), len(M), len(T), len(ITER), len(S)))

start = time.time()
for itt in ITER:
    dsp_sum = 0
    for ss in S:
        DSP_Obj_Value, chi_plus_sol, chi_minus_sol, GAMMA_plus_sol, GAMMA_minus_sol = Dual_Sub_Problem(ss,
                     I, J, K, T, M,
                     ca2, ca3, ca4, c2, c3, c4, c5, c6,
                     v1, v2, alpha, dem, x_star, landa_star, ca4_upper, ca4_lower, dem_upper, dem_lower)
        
        dsp_sum += DSP_Obj_Value
        for i in I:
            for t in T:
                for k in K:
                    chi_plus_star[i,k,t,itt,ss] = chi_plus_sol[i,k,t]
                    chi_minus_star[i,k,t,itt,ss] = chi_minus_sol[i,k,t]
                for m in M:
                    GAMMA_plus_star[i,m,t,itt,ss] = GAMMA_plus_sol[i,m,t]
                    GAMMA_minus_star[i,m,t,itt,ss] = GAMMA_minus_sol[i,m,t]

    ub = (quicksum(f1[i,j,n]*zeta_star[i,j,n] for i in I for j in J for n in N) + quicksum(c1[i,j,k]*x_star[i,j,k] for i in I for j in J for k in K) + (1/len(S))*dsp_sum) + landa_star*theta
    UB = ub.getValue()
    print(f"UB for iteration {itt}: {UB}")
    rmp_value, x_sol, zeta_sol, landa_sol = Restricted_Master_Problem(itt,
                              I, J, K, N, S, T, M,
                              f1, c1, g, ca1, ca2, ca3, ca4, 
                              u, BigM, theta, ca4_upper, ca4_lower, dem_upper, dem_lower,
                              chi_plus_star, chi_minus_star, GAMMA_plus_star, GAMMA_minus_star)

    LB = rmp_value

    landa_star = landa_sol
    for i in I:
        for j in J:
            for n in N:
                zeta_star[i,j,n] = zeta_sol[i,j,n]
            for k in K:
                x_star[i,j,k] = x_sol[i,j,k]

    
    Gap = UB - LB
    time_now = time.time()
    
    print(f"LB for iteration {itt}: {LB}")
    print(f"Gap for iteration {itt}: {Gap}")
    if (Gap < 0.01) or (time_now - start > 3600):
        break
    print(40*'-')
end = time.time()
print(f"Time taken: {end - start}")