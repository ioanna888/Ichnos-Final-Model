"""Παράγει τον μάρτυρα χωρίς ανάδραση + το πλήρες CSV, με roadrunner για το
τρέχον κύκλωμα και ισοδύναμη ODE υλοποίηση για τον μάρτυρα (δομική αλλαγή
που δεν γίνεται με απλή αλλαγή παραμέτρου)."""
import csv, numpy as np, roadrunner
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

b,u_w,delta_w,a_TetR,k_deg_TetR,K_R,n,P_min,mu = 4,1,0.12,50,0.12,0.44,4,0.0056,0.35
a_Rep,m,m1,m2,k_deg_Rep,E,f,eps = 50,7.39,2.46,1.37,0.01,0.173,1.0,1e-9
k_deg_TIP=1.0
def beta(S,EC50,n_s,bb=5,bm=25): return bb+(bm-bb)*(S/EC50)**n_s/(1+(S/EC50)**n_s)

def run(S,EC50,n_s,feedback=True,P_fixed=None):
    B=beta(S,EC50,n_s)
    def rhs(t,y):
        TIP,TetR,C,Rdg,Rg,Rdr,Ri,Rr=y
        P=P_min+(1-P_min)/(1+(max(TetR,0)/K_R)**n)
        Pt = P if feedback else P_fixed
        bind=b*TIP*TetR-u_w*C
        return [B-bind-(k_deg_TIP+mu)*TIP, a_TetR*Pt-bind-(mu+k_deg_TetR)*TetR, bind-(delta_w+mu)*C,
                a_Rep*P-m*Rdg-(k_deg_Rep+mu)*Rdg, m*Rdg-(k_deg_Rep+mu)*Rg,
                a_Rep*P-m1*Rdr-(k_deg_Rep+mu)*Rdr, m1*Rdr-m2*Ri-(k_deg_Rep+mu)*Ri, m2*Ri-(k_deg_Rep+mu)*Rr]
    s=solve_ivp(rhs,[0,4000],[0]*8,rtol=1e-11,atol=1e-13,method='LSODA')
    TIP,TetR,C,Rdg,Rg,Rdr,Ri,Rr=s.y[:,-1]
    bf=Rr/(Rdr+Ri+Rr+eps)
    return dict(TIP=TIP,TetR=TetR,C=C,OG=Rg*(1-bf*E))

def hill4(S,ba,am,K,nn): return ba+am*S**nn/(K**nn+S**nn)
def met(S,y):
    p,_=curve_fit(hill4,S,y,p0=[y.min(),y.max()-y.min(),np.median(S),1.5],maxfev=400000)
    A=np.vstack([S,np.ones_like(S)]).T; c,*_=np.linalg.lstsq(A,y,rcond=None)
    return abs(p[3]), 1-(y-A@c).var()/y.var(), y.max()/y.min()

GRID={"ox":(271,1.7,np.array([10.,20,40,80,150,271,400,600])),
      "er":(2200,2.0,np.array([130.,260,500,900,1500,2200,2800,3300]))}
SID={"ox":"S_ox","er":"S_er"}

rows=[]
for v,(EC50,n_s,Sg) in GRID.items():
    # τρέχον κύκλωμα — από το ΠΡΑΓΜΑΤΙΚΟ merged (roadrunner)
    yr=[]
    for S in Sg:
        r=roadrunner.RoadRunner(f"exportsbml/merged_{v}.sbml"); r.resetToOrigin()
        r[SID[v]]=float(S); yr.append(r.simulate(0,400,4000,["time","Observed_Green"])[-1,1])
    yr=np.array(yr); n1,r1,f1=met(Sg,yr)
    # μάρτυρας χωρίς ανάδραση
    Pf=run(Sg[0],EC50,n_s,True)['P'] if False else (P_min+(1-P_min)/(1+(run(Sg[0],EC50,n_s,True)['TetR']/K_R)**n))
    y2=np.array([run(S,EC50,n_s,False,Pf)['OG'] for S in Sg]); n2,r2,f2=met(Sg,y2)
    # είσοδος
    yb=np.array([beta(S,EC50,n_s) for S in Sg]); nb,rb,fb=met(Sg,yb)
    for lab,(ne,rr,fc) in [("ΕΙΣΟΔΟΣ beta(S)",(nb,rb,fb)),("ΧΩΡΙΣ ανάδραση (μάρτυρας)",(n2,r2,f2)),
                           ("ΜΕ ανάδραση (ICHNOS)",(n1,r1,f1))]:
        rows.append(dict(variant=v,circuit=lab,n_eff=round(ne,4),R2_linear=round(rr,4),fold_change=round(fc,4)))
    print(f"{v}: χωρίς={n2:.4f}/{r2:.4f}/{f2:.4f}   με={n1:.4f}/{r1:.4f}/{f1:.4f}   "
          f"Δn_eff={100*(1-n1/n2):.1f}%")

with open("control_results.csv","w",newline="",encoding="utf-8") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("γράφτηκε control_results.csv")
