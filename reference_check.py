# reference_check.py - independent NumPy implementation of the same discrete equations as the Julia solvers.
# Runs every validation check and produces validation_results.png.   Usage: python reference_check.py
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
res={}
def jac1d(N,u,L,rho,g,scheme,eps=1e-7):
    h=L/N; x=np.linspace(h/2,L-h/2,N); F=rho*u; D=g/h
    if scheme=="central": aE,aW,a0,aL=D-F/2,D+F/2,2*D+F,2*D-F
    else: aE,aW,a0,aL=D,D+F,2*D+F,2*D
    phi=np.zeros(N); r=1; it=0; hist=[]
    while r>eps:
        n=phi.copy()
        n[1:-1]=(aE*phi[2:]+aW*phi[:-2])/(aE+aW)
        n[0]=(aE*phi[1]+a0*1.0)/(aE+a0)
        n[-1]=(aW*phi[-2]+aL*0.0)/(aW+aL)
        r=np.abs(n-phi).sum(); phi=n; it+=1; hist.append(r)
    A=np.zeros((N,N)); b=np.zeros(N)
    for i in range(N):
        if i==0: A[i,i]=aE+a0; A[i,1]=-aE; b[i]=a0
        elif i==N-1: A[i,i]=aW+aL; A[i,i-1]=-aW
        else: A[i,i]=aE+aW; A[i,i+1]=-aE; A[i,i-1]=-aW
    ref=np.linalg.solve(A,b)
    Pe=F*L/g; xf=np.linspace(0,L,400); ex=1-(np.exp(Pe*xf)-1)/(np.exp(Pe)-1)
    exn=1-(np.exp(Pe*x)-1)/(np.exp(Pe)-1)
    return x,phi,ref,xf,ex,exn,it,hist,F/D
xc,pc,rc,xf,ex,exc,itc,hc,pec=jac1d(20,2.5,1,1,0.1,"central")
xu,pu,ru,_,_,exu,itu,hu,peu=jac1d(5,2.5,1,1,0.1,"upwind")
res["1D central (N=20) vs direct solve"]=np.abs(pc-rc).max()
res["1D central (N=20) vs analytic"]=np.abs(pc-exc).max()
res["1D upwind (N=5) vs direct solve"]=np.abs(pu-ru).max()
res["1D upwind (N=5) vs analytic"]=np.abs(pu-exu).max()
res["1D upwind bounded in [0,1]"]=bool(((pu>=0)&(pu<=1)).all())
ref_err={}
for Nn in [20,40,80,160]:
    o=jac1d(Nn,2.5,1,1,0.1,'central',eps=1e-10); ref_err[Nn]=np.abs(o[2]-o[5]).max()
print('refinement',ref_err)
# conduction
N=20;x=np.linspace(0,1,N);T=np.zeros(N);T[-1]=1;r=1;itk=0
while r>1e-7:
    n=T.copy(); n[1:-1]=0.5*(T[2:]+T[:-2]); r=np.abs(n-T).sum(); T=n; itk+=1
res["1D conduction vs linear exact"]=np.abs(T-x).max()
# 2D (same coefficients as TwoDim_Conv_Diff.jl)
N=25;L=1;h=L/N;g=1e-4;F=5.0;Fy=1.0;D=g/h;pL,pR,pT,pB=1,0,1,0
phi=np.zeros((N,N));err=1;it2=0
def coeffs(i,j):
    aW=D+F if j>0 else 0; aE=D if j<N-1 else 0; aN=D if i>0 else 0; aS=D+Fy if i<N-1 else 0
    Su=0;aP=aW+aE+aN+aS
    if j==0: Su+=(2*D+F)*pL; aP+=2*D+F
    if j==N-1: Su+=2*D*pR; aP+=2*D
    if i==0: Su+=2*D*pT; aP+=2*D
    if i==N-1: Su+=(2*D+Fy)*pB; aP+=2*D+Fy
    return aW,aE,aN,aS,Su,aP
C=[[coeffs(i,j) for j in range(N)] for i in range(N)]
while err>1e-6:
    n=np.empty_like(phi)
    for i in range(N):
        for j in range(N):
            aW,aE,aN,aS,Su,aP=C[i][j]
            s=Su
            if j>0:s+=aW*phi[i,j-1]
            if j<N-1:s+=aE*phi[i,j+1]
            if i>0:s+=aN*phi[i-1,j]
            if i<N-1:s+=aS*phi[i+1,j]
            n[i,j]=s/aP
    err=np.abs(n-phi).sum(); phi=n; it2+=1
A=np.zeros((N*N,N*N));b=np.zeros(N*N);idx=lambda i,j:i*N+j
for i in range(N):
    for j in range(N):
        aW,aE,aN,aS,Su,aP=C[i][j];k=idx(i,j);A[k,k]=aP;b[k]=Su
        if j>0:A[k,idx(i,j-1)]=-aW
        if j<N-1:A[k,idx(i,j+1)]=-aE
        if i>0:A[k,idx(i-1,j)]=-aN
        if i<N-1:A[k,idx(i+1,j)]=-aS
ref2=np.linalg.solve(A,b).reshape(N,N)
res["2D conv-diff (25x25) vs direct solve"]=np.abs(phi-ref2).max()
res["2D bounded in [0,1]"]=bool(((phi>-1e-6)&(phi<1+1e-6)).all())
# FD demo
f=lambda x:-4*x**3+7*x**2-3*x+9; x0,hh=1.5,0.1; ex_d=-12*x0**2+14*x0-3
fd=(f(x0+hh)-f(x0))/hh; bd=(f(x0)-f(x0-hh))/hh; cd=(f(x0+hh)-f(x0-hh))/(2*hh)
res["Central diff error = -4h^2 (truncation)"]=abs((cd-ex_d)-(-4*hh**2))
for k,v in res.items(): print(f"{k:42s} {v}")
print("iters",itc,itu,itk,it2,"cellPe",pec,peu, "fd/bd/cd/exact",fd,bd,cd,ex_d)
# figure
fig,ax=plt.subplots(2,3,figsize=(15,8.5))
a=ax[0,0];a.plot(xf,ex,'k-',label='Analytic');a.plot(xc,pc,'s',mfc='none',label='Central, N=20');a.plot(xu,pu,'o',label='Upwind, N=5');a.set_title('1D advection–diffusion (Pe = 25)');a.set_xlabel('x');a.set_ylabel('φ');a.legend()
a=ax[0,1];a.semilogy(hc,label=f'Central N=20 ({itc} it)');a.semilogy(hu,label=f'Upwind N=5 ({itu} it)');a.set_title('Jacobi convergence');a.set_xlabel('iteration');a.set_ylabel('residual');a.legend()
a=ax[0,2];Ns=np.array(list(ref_err));Es=np.array(list(ref_err.values()));a.loglog(Ns,Es,'o-',label='Central scheme error');a.loglog(Ns,Es[0]*(Ns[0]/Ns)**2,'k--',label='slope 2 (2nd order)');a.set_title('Grid refinement vs analytic');a.set_xlabel('N cells');a.set_ylabel('max |φ − φ_exact|');a.legend()
a=ax[1,0];yc=L-np.linspace(h/2,L-h/2,N);xx=np.linspace(h/2,L-h/2,N);cs=a.contourf(xx,yc,np.clip(phi,0,1),levels=np.linspace(0,1,41),cmap='viridis');fig.colorbar(cs,ax=a);a.set_title(f'2D convection–diffusion, u=5, v=1 ({it2} it)');a.set_xlabel('x');a.set_ylabel('y');a.set_aspect('equal')
a=ax[1,1];cs=a.contourf(xx,yc,np.abs(phi-ref2),levels=40,cmap='magma');fig.colorbar(cs,ax=a);a.set_title('2D |iterative − direct|');a.set_aspect('equal')
a=ax[1,2];a.axis('off')
rows=[(k,(f"{v:.1e}" if isinstance(v,float) else ("PASS" if v else "FAIL"))) for k,v in res.items()]
t=a.table(cellText=rows,colLabels=['Check','Max error / result'],loc='center',colWidths=[0.75,0.3]);t.auto_set_font_size(False);t.set_fontsize(8.5);t.scale(1,1.5);a.set_title('Validation summary')
fig.suptitle('CFD solver validation – Computational-Fluid-Dynamics-Julia',fontsize=13);fig.tight_layout()
fig.savefig('validation_results.png',dpi=130)
