from gurobipy import *
import numpy as np
import random
import time

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
ITER = range(100)

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


def Dual_Sub_Problem(sample,
                     I, J, K, T, M,
                     ca2, ca3, ca4, c2, c3, c4, c5, c6,
                     v1, v2, alpha, dem, x_star):
    DSP = Model("Dual Sub Problem")
    DSP.setParam('OutputFlag', 0)

    gamma = DSP.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="gamma")
    mu = DSP.addVars(I, K, T, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="mu")
    eta = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="eta")
    nu = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="nu")
    rho = DSP.addVars(I, M, T, lb=0, vtype=GRB.CONTINUOUS, name="rho")
    pi = DSP.addVars(I, K, T, lb=0, vtype=GRB.CONTINUOUS, name="pi")

    obj_fn = (-1*quicksum(ca4[i,k,sample]*nu[i,k,t] for i in I for k in K for t in T) 
            - quicksum((ca2[i,k]-sum(x_star[i,j,k] for j in J))*eta[i,k,0] for i in I for k in K) 
            - quicksum(ca2[i,k]*eta[i,k,t] for i in I for k in K for t in T if t != 0) 
            - quicksum(x_star[i,j,k]*gamma[i,k,0] for i in I for j in J for k in K) 
            - quicksum(ca3[i,k]*pi[i,k,t] for i in I for k in K for t in T if t != len(T) - 1) 
            - quicksum(dem[i,m,t,sample]*rho[i,m,t] for i in I for m in M for t in T))
    
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
    
    DSP.optimize()


    if DSP.status == GRB.OPTIMAL:

        DSP_Obj_Value = DSP.ObjVal
        gamma_sol = {(i,k,t): gamma[i,k,t].X for i in I for k in K for t in T}
        mu_sol = {(i,k,t): mu[i,k,t].X for i in I for k in K for t in T}
        eta_sol = {(i,k,t): eta[i,k,t].X for i in I for k in K for t in T}
        nu_sol = {(i,k,t): nu[i,k,t].X for i in I for k in K for t in T}
        pi_sol = {(i,k,t): pi[i,k,t].X for i in I for k in K for t in T}
        rho_sol = {(i,m,t): rho[i,m,t].X for i in I for m in M for t in T}

        return DSP_Obj_Value, gamma_sol, mu_sol, eta_sol, nu_sol, pi_sol, rho_sol
    else:
        print("❌ DSP ERROR ❌")  
        return None, None, None, None, None, None, None
    



  
def Restricted_Master_Problem(iter,
                              I, J, K, N, S, T, M,
                              f1, c1, g, ca1, ca2, ca3, ca4, 
                              u, BigM,
                              nu_star, eta_star, gamma_star, pi_star, rho_star):
    IT = range(iter+1)
    RMP = Model("Restricted Master Problem")
    RMP.setParam('OutputFlag', 0)

    delta = RMP.addVars(J, vtype=GRB.BINARY, name="delta")
    zeta = RMP.addVars(I, J, N, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="zeta")
    x = RMP.addVars(I, J, K, lb=0, vtype=GRB.CONTINUOUS, name="x")
    phi_cut = RMP.addVars(S, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="phi_cut")

    obj_fn = (
    quicksum(f1[i,j,n]*zeta[i,j,n] for i in I for j in J for n in N)
    + quicksum(c1[i,j,k]*x[i,j,k] for i in I for j in J for k in K)
    + (1/len(S))*quicksum(phi_cut[s] for s in S))

    RMP.setObjective(obj_fn, GRB.MINIMIZE)

    for i in I:
        for j in J:
            RMP.addSOS(GRB.SOS_TYPE2, [zeta[i, j, n] for n in N])

    for i in I:
        for j in J:
            # --- Constraint 1 ---
            RMP.addConstr(quicksum(x[i,j,k] for k in K) <= ca1[i,j])

            # --- Constraint 3 ---
            RMP.addConstr(quicksum(x[i,j,k] for k in K) == quicksum(g[i,j,n]*zeta[i,j,n] for n in N))

            # --- Constraint 5 ---
            for k in K:
                RMP.addConstr(x[i,j,k] <= BigM*delta[j])
            
            # --- Constraint 6 ---
            RMP.addConstr(quicksum(zeta[i,j,n] for n in N) == 1)

    # --- Constraint 2 ---
    for i in I:
        for k in K:
            RMP.addConstr(quicksum(x[i,j,k] for j in J) <= ca2[i,k])

    # --- Constraint 4 ---
    for j in J:
        for jp in J:
            RMP.addConstr(delta[j] + delta[jp] <= 1 + u[j,jp])

    # --- Optimality Cut ---
    for iter in IT:
        for sampel in S:
            lhs = phi_cut[sampel]
            rhs = (-1*quicksum(ca4[i,k,sampel]*nu_star[i,k,t,iter,sampel] for i in I for k in K for t in T) 
            - quicksum((ca2[i,k]-sum(x[i,j,k] for j in J))*eta_star[i,k,0,iter,sampel] for i in I for k in K) 
            - quicksum(ca2[i,k]*eta_star[i,k,t,iter,sampel] for i in I for k in K for t in T if t != 0) 
            - quicksum(x[i,j,k]*gamma_star[i,k,0,iter,sampel] for i in I for j in J for k in K) 
            - quicksum(ca3[i,k]*pi_star[i,k,t,iter,sampel] for i in I for k in K for t in T if t != len(T) - 1) 
            - quicksum(dem[i,m,t,sampel]*rho_star[i,m,t,iter,sampel] for i in I for m in M for t in T))
            RMP.addConstr(lhs >= rhs)
    
    RMP.optimize()

    if RMP.status == GRB.OPTIMAL:
        RMP_Obj_Value = RMP.ObjVal

        delta_sol = {j: delta[j].X for j in J}
        x_sol = {(i,j,k): x[i,j,k].X for i in I for j in J for k in K}
        zeta_sol = {(i,j,n): zeta[i,j,n].X for i in I for j in J for n in N}
        
        return RMP_Obj_Value, delta_sol, x_sol, zeta_sol

    else:
        print("❌ RMP ERROR ❌")  
        return None, None, None, None

delta_star = np.zeros((len(J)))
zeta_star = np.zeros((len(I), len(J), len(N)))
x_star = np.zeros((len(I), len(J), len(K)))

num_random_elements = int(0.1 * len(I) * len(J) * len(K))

indices = set()
while len(indices) < num_random_elements:
    i = np.random.randint(0, len(I))
    j = np.random.randint(0, len(J))
    k = np.random.randint(0, len(K))
    indices.add((i, j, k))

for i, j, k in indices:
    x_star[i, j, k] = np.random.uniform(5, 15)

gamma_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
mu_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
eta_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
nu_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))
rho_star = np.zeros((len(I), len(M), len(T), len(ITER), len(S)))
pi_star = np.zeros((len(I), len(K), len(T), len(ITER), len(S)))


for itt in ITER:
    dsp_sum = 0
    for ss in S:
        dsp_value, gamma_sol, mu_sol, eta_sol, nu_sol, pi_sol, rho_sol = Dual_Sub_Problem(ss,
                        I, J, K, T, M,
                        ca2, ca3, ca4, c2, c3, c4, c5, c6,
                        v1, v2, alpha, dem, x_star)
        for i in I:
            for t in T:
                for k in K:
                    gamma_star[i,k,t,itt,ss] = gamma_sol[i,k,t]
                    mu_star[i,k,t,itt,ss] = mu_sol[i,k,t]
                    eta_star[i,k,t,itt,ss] = eta_sol[i,k,t]
                    nu_star[i,k,t,itt,ss] = nu_sol[i,k,t]
                    pi_star[i,k,t,itt,ss] = pi_sol[i,k,t]
                for m in M:
                    rho_star[i,m,t,itt,ss] = rho_sol[i,m,t]

        dsp_sum += dsp_value

    ub = quicksum(f1[i,j,n]*zeta_star[i,j,n] for i in I for j in J for n in N) + quicksum(c1[i,j,k]*x_star[i,j,k] for i in I for j in J for k in K) + (1/len(S))*dsp_sum
    UB = ub.getValue()
    print(f"UB for iteration {itt}: {UB}")

    rmp_value, delta_sol, x_sol, zeta_sol = Restricted_Master_Problem(itt,
                                I, J, K, N, S, T, M,
                                f1, c1, g, ca1, ca2, ca3, ca4, 
                                u, BigM,
                                nu_star, eta_star, gamma_star, pi_star, rho_star)
    LB = rmp_value
    for i in I:
        for j in J:
            for n in N:
                zeta_star[i,j,n] = zeta_sol[i,j,n]
            for k in K:
                x_star[i,j,k] = x_sol[i,j,k]
    
    print(f"LB for iteration {itt}: {LB}")
    Gap = UB - LB
    print(f"Gap for iteration {itt}: {Gap}")

    if (Gap < 0.01):
        break
    print(40*'-')