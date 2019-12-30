# -*- coding: utf-8 -*-
# rl = 'import Post\nreload( Post)\nfrom  Post import *'
"""
Created on Fri Dec 27 16:10:43 2019

Sage implementation of the algorithms in
"Complete algorithms for algebraic strongest postconditions and
weakest preconditions in  polynomial ODEs"
by Michele Boreale, December 2019.

Author: Michele Boreale, University of Florence, Italy.

IMPORTANT: this code must be executed in Python *under SageMath*. Tested with SageMath 8.1 under Windows 10.
Import the script from the Python console with:

from Post import *

1. MAIN PROCEDURES
See also the paper "Complete algorithms for algebraic strongest postconditions and
weakest preconditions in  polynomial ODEs", by Michele Boreale, 2019.


1.1 out = post(pt,F,P,Xlist,Plist)

where

INPUT arguments are
    - pt: a polynomial template (possibly generated with genpt, see 2 below)
    - F: a list of polynomials in the variables in Xlist, representing a vector field
    - P: a list of polynomials or an ideal, representing the precondition psi
    - Xlist: list of variables
    - Plist: list of parameters in pt

OUTPUT out can be one of
    - False: pt is not an invariant (this may arise only in case pt is actually a polynomial)
    - (ptout,J), where:
        * ptout is the result polynomial template (the most general instance of pt all of which instances are valid polynomial invariants)
        * J is the smallest invariant ideal containing all the valid instances of pt; var(J) is the largest algebraic (invariant) precondition implying these laws
    - (ptout,G,"Warning..."), like above, but with a warning that the maximum n. of iterations has been reached (result might not be correct).


1.2 out = pre(F,P,Xlist)

where

INPUT arguments are
    - F: a list of polynomials in the variables in Xlist, representing a vector field
    - P: a list of polynomials, representing the postcondition phi
    - Xlist: a list of variables

OUTPUT out is an ideal I such that var(I) is the weakest precondition for phi.

2. ANOTHER USEFUL PROCEDURE

1) pt,plist = genpt(Xlist,deg)
INPUT
    - Xlist: a list of variables
    - deg: a nonnegative integer

OUTPUT
    - pt: a complete polynomial template of degree deg built out of variables in Xlist
    - plist: the corresponding list of parameters (one for each monomial)



3. EXAMPLES IN THE PAPER
NB: the necessary variables are already defined in the code.

1) Collision avoidance.
Try the following.

pt,pl=genpt(xColl,2)                                            # generates a template and a parameter list
P=[x1-x10,y1-y10,x2-x20,y2-y20,d1-d10,d2-d20,e1-e10,e2-e20]     # 'generic' precondition psi
outcoll = post(pt,COLL,P,xColl,pl) 
J=outcoll[1]              # smallest invariant ideal containing all valid instances of pt
J.groebner_basis()        # a 12 elements basis
J.ring()[pl](outcoll[0])  # prints a readable version of the result template


2) Airplane vertical motion.
Try the following.

pt,pl=genpt(xAir+[q*u,q*w],2)                                   # generates a template and a parameter list
P=[theta,sinth,costh-1,u-u0,w-w0,x-x0,z-z0,q-q0]                # precondition psi0
outair = post(pt,AIR,P,xAir,pl) 
J=outair[1]              # smallest invariant ideal implying all valid instances of pt
J.groebner_basis()       # a 15 elements basis
J.ring()[pl](outair[0])  # prints a readable version of the result template



3) Kepler laws.
Try the following.

a1=var('a1')
pt1=a1*pell                                                      # generate simple template for elliptic orbits
P0=[a-1,e, GM-1,r-1,vr,u-1,theta,omega-1]                        # precondition psi0 implying circular orbit of radius 1
outkepl1 = post(pt1,KEPL,P0,xKepl,[a1])
J = outkepl1[1]                                                  # J = smallest invariant ideal implying pell=0 (1st Kepler law)
J0 = (J.ring()).ideal([dA-r**2*omega/2, r-a*(1-e), theta, vr, u*r-1, costh-1, sinth, a*b**2-1, e-c**2, (1-e)*d**2-1, GM*f**2-1])    # J0 encodes consistency and positivity conditions
JR = (J+J0).radical()                                            # V(JR) gives largest *physically meaningful* algebraic invariant variety implying pell=0
I = (J.ring()).ideal([r**2*omega**2-GM*u*(e+1), dA-r**2*omega/2, r-a*(1-e), theta, vr, u*r-1, costh-1, sinth, a*b**2-1, e-c**2, (1-e)*d**2-1, GM*f**2-1]) 
                                                                 # I is the ideal defined by P in (18) in the paper
I==JR                                                            # This gives True. Hence psi=V(I)=V(JR).
pt2,pl2 = genpt([GM,a,e,r,u,dA],4)                               # generates a template and a parameter list over a subset of the variables
outkepl23 = post(pt2,KEPL,I,xKepl,pl2)                           # outkepl23[0] is the result template, it contains 2nd and 3rd Kepler law as instances
J.ring()[pl2](outkepl23[0])    # prints a readable version of the result template; look for summands with 2nd Kepler and 3rd Kepler laws
outkepl23[1].groebner_basis() # 2 elements Groebner basis of smallest invariant ideal implying 2nd and 3rd K. laws.


4) Running example.
pt,pl = genpt(xRun,2)
outrunpost = post(pt,RUN,[x-y],xRun,pl)
outrunpre  = pre(RUN,[x**2-x*y],xRun)


5) 3D Lotka-Volterra
x0,y0,z0,x,y,z = var('x0,y0,z0,x,y,z')
xLV=[x0,y0,z0,x,y,z]
LV = [0,0,0,x*y-x*z,y*z-y*x,z*x-z*y]   # 3D Lotka-Volterra vector field
pt,pl=genpt(xLV,3)
PLV=[x-x0,y-y0,z-z0]
qt,JLV=post(pt,LV,PLV,xLV,pl)
JLV.groebner_basis()

6) Coupled spring-mass
kom,L,x10,x20,x1,x2,v1,v2 = var('kom,L,x10,x20,x1,x2,v1,v2')
xSpring = [kom,L,x10,x20,x1,x2,v1,v2]
SPRING  = [0,0,0,0,v1,v2, kom*(x2-2*x1),-kom*(x2-x1-L)]
pt,pl=genpt(xSpring,3)
PSpring=[x1-x10,x2-x20,v1,v2]
qt,JSpring=post(pt,SPRING,PSpring,xSpring,pl)
JSpring


"""

from sage.all import *
from sympy  import itermonomials as itm
import time



####   DEFINITIONS OF EXAMPLES
#1. Collision Avoidance 

omega1,omega2,x10,x20,y10,y20,d10,d20,e10,e20,x1,x2, y1,y2,d1,d2,e1,e2 = var('omega1,omega2,x10,x20,y10,y20,d10,d20,e10,e20,   x1,x2, y1,y2,    d1,  d2,     e1,            e2')

xColl = [omega1,omega2,x10,x20,y10,y20,d10,d20,e10,e20,x1,x2,y1,y2,d1,d2,e1,e2]

COLL =  [0,0,0,0,0,0,0,0,0,0,d1,d2,e1,e2,-omega1*d2,omega1*d1,-omega2*e2,omega2*e1]


# 2. Airplane vertical motion

u,w,x,z,q,theta,costh,sinth,g,Xom,Zom,MoIyy,u0,w0,x0,z0,q0 = var('u,w,x,z,q,theta,costh,sinth,g,Xom,Zom,MoIyy,u0,w0,x0,z0,q0')

xAir = [u,w,x,z,q,theta,costh,sinth,g,Xom,Zom,MoIyy,u0,w0,x0,z0,q0]

AIR = [Xom-g*sinth-q*w,Zom+g*costh+q*u,costh*u+sinth*w,-sinth*u+costh*w,MoIyy,q,-sinth*q,costh*q,0,0,0,0,0,0,0,0,0]


# 3. Kepler laws
a,b,c,d,e,f,GM,r,theta,vr,omega,u,costh,sinth,dA=var('a,b,c,d,e,f,GM,r,theta,vr,omega,u,costh,sinth,dA')

xKepl = [a,b,c,d,e,f,GM,r,theta,vr,omega,u,costh,sinth,dA]

KEPL = [0,0,0,0,0,0,0,vr,omega,-GM*u**2+r*omega**2,-2*vr*omega*u,-u**2*vr,-sinth*omega,costh*omega,-omega*r**2*u*vr + omega*r*vr]

pell = r+r*e*costh-a*(1-e**2)

assumptpolalt = [dA-r**2*omega/2, r-a*(1-e), theta, vr, u*r-1, costh-1, sinth, a*b**2-1, e-c**2, (1-e)*d**2-1, GM*f**2-1]
omegainalt =  r**2*omega**2-GM*u*(e+1)             
                

K1=pell
K2= -omega*r**2*u*vr+omega*r*vr 
K3=-(e**2 - 1)*GM*a - 4*dA**2


# 4. Running example

x,y = var('x,y')
xRun= [x,y]
RUN = [y**2,x*y]


##################################################################

maxiter = 60
expandflag = True
algebraiclosure = False
verbose = False

def genpt(Xlist,deg):
    monlist=list(itm(Xlist,deg))
    l=len(monlist)
    parlist = list(var('a%d' % j) for j in range(l))
    prod=((Matrix(monlist))*Matrix(parlist).T)[0,0]
    return prod, parlist
    
def Lie(pt,F,Xlist):
    global VF # vector field as a column vector
    return expansion((Matrix(pt.gradient(Xlist))*VF)[0,0])  # Lie der.    
    
def stringify(L):
    s = "".join(str(x)+',' for x in L)
    return  s[0:len(s)-1]

def linind(list,p):     # stub: should check linear independence of p from res
    return not(p==0)

def instantiate(ptlist):
    global zeropar,R
    res = []
    for pt in ptlist:
        if zeropar.keys()==[]:
            if linind(res,pt):
                res.append(pt)
        else:
            for a in zeropar.keys():
                tmp=zeropar.copy()
                tmp.update({a:1})
                p=R(pt.subs(tmp))
                if linind(res,p):
                    res.append(p)
    return res # returns list of polynomials in R

def included(newinst,H):
    for p in newinst:
        if (not(p in H)):
            return False
    return True

def seqsub(sigma,pt):   # sequentially applies a list og substitutions in sigma to template pt
    for s in sigma:
        pt=pt(s)
    return pt

def expansion(p):
    if expandflag:
        return expand(p)
    else:
        return p

def initglobal(pt,F,P,Xlist,Plist):
    global start_time, VF, K, R, RU, Rpar, I0, zeropar, freepar, GS, uptodateGS
    singular.quit()
    if algebraiclosure:
        K=QQbar
    else:
        K=QQ
    start_time = time.time()
    VF = Matrix(F).T  # vector field as a column vector
    xstr = stringify(Xlist)
    pstr = stringify(Plist)
    pxstr= stringify(Plist+Xlist)
    R = K[xstr] # polynomial ring
    RU = K[pxstr] # universal polynomial ring
    if pstr == '':
        Rpar = R
    else:
        Rpar=(K[pstr])[xstr] # polynomial ring of templates
    if type(P)==type([]):
        I0=RU.ideal(P).radical()               # radical of <x0genlist>, again in RU
    else:
        I0=P.radical()
    freepar = set(pt.variables())-set(Xlist) 
    zeropar = {a:0 for a in freepar}
    uptodateGS = False
    GS = []
    
                
def reporterror(j,gamma):
    global start_time                
    print("Linear constraint for V_"+str(j)+": "+str(gamma))
    print("Cannot find instance of given template that is an invariant.")
    print('m='+str(j))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))


def solver(rt):
    global Rpar, freepar,zeropar, uptodateGS
    gamma=[]
    coeffs = [SR(lin) for lin in Rpar(SR(rt)).coefficients()] # list of linear expressions extracted as coefficients of rt                                                   
    for lin in coeffs:      # solve the list of linear constraints, one by one
        lin=expansion(seqsub(gamma,lin))  # apply previously accumulated subst. to current constraint lin
        if lin.variables()==():
            if lin==SR(0):
                C=True
            else:
                gamma=gamma+[{lin:0}]
                return False, gamma
        else:
            C=solve(lin,lin.variables()[0],solution_dict=True)      # solve one linear constraint lin
            s = C[0]                                                # obtain subst. s from solution
            gamma=gamma+[s]
    freepar = freepar-set([ s.keys()[0] for s in gamma])
    zeropar = { a:0 for a in freepar }
    if gamma!=[]:    
        uptodateGS=False                                         # append s to list of substitutions
    return True, gamma

        
def printiter(j,pt,rt):
    print("j="+str(j))
    if verbose:
        print("  pt"+str(j)+"= ",pt)
        print("  r"+str(j)+"= ",rt)
    else:
        print("  pt"+str(j)+" and r"+str(j)+" computed.")

def reportstable(j):
    global start_time
    print("  Equality holds, both chains stabilized")
    print('Stabilization at m='+str(j-1))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))      
    
    
def post(pt,F,P,Xlist,Plist):
    global start_time, VF, K, R, RU, Rpar, I0, zeropar, freepar, GS, uptodateGS
    initglobal(pt,F,P,Xlist,Plist)                                         
    
    S=[]       
    for j in range(maxiter):        
        rt=RU(I0.reduce(RU(pt)))     # rt = reduce pt modulo I0
        flag,gamma=solver(rt)
        printiter(j,pt,rt)
        if flag==False:             # report error if rt=0 is unsatisfiable (may occur only if argument pt is not a proper template)
            reporterror(j,gamma)
            return False        
        print("  Linear constraints for V"+str(j)+": "+str(gamma))
        if gamma==[]:               # Vj+1 = Vj
            print("  Vector space equality detected: V"+str(j)+" = V"+str(j-1))
            print("  Checking ideal equality J"+str(j)+" = J"+str(j-1)+"...")
            if not(uptodateGS):
                GS = instantiate(S)
                uptodateGS = True
            J=R.ideal(GS)
            ptinst=instantiate([pt])
            if included(ptinst,J):    # Ideals stabilized => termination
                reportstable(j)              
                return S[0],J
            else:
                print("  Equality does not hold, chains not yet stabilized")
                GS = GS + ptinst
        else:                     # Ideal chain not stabilized
            pt = SR(expansion(seqsub(gamma,pt)))       # propagate constraints gamma on pt
            S = [SR((seqsub(gamma,qt))) for qt in S]   # propagate constraints gamma on S
        if not(pt==SR(0)):
            S.append(pt)
        pt=Lie(pt,F,Xlist)       # compute Lie-derivative of pt
    print('  m='+str(j))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))    
    return S[0],GS,'WARNING: maximum number of iterations reached, stabilization not reached.'
        

############## PRE algorithm ########################

MAXITER2 = 50        

def pre(F,P,Xlist): 
    singular.quit()
    singular.lib('poly.lib')
    singular.lib('primdec.lib')
    start_time = time.time()
    Xstring = stringify(Xlist)
    polylist=P
    singular.eval('ring  R = 0, ('+Xstring+'),lp')
    singular.eval('ideal J ='+stringify(polylist))
    singular.eval('ideal JR = radical(J)')
    
    j=0
    while j<MAXITER2:
        print(j)
        #print(polylist)
        newptlist=[]
        singular.eval('ideal GR=groebner(JR);' )
        for pt in P:
            newpt=Lie(pt,F,Xlist)
            newptst=str(newpt)
            singular.eval('poly r = reduce('+newptst+',GR)')
            r=singular.poly('r').sage()
            #print(r)
            if SR(r)!=0:
                newptlist.append(SR(newpt))            
        if newptlist==[]:
            print("Chain stabilized")
            print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
            I = singular.ideal('GR')
            return I.groebner()        

        polylist=polylist+newptlist
        singular.eval('ideal JR=' +stringify(polylist)+';')
        P=newptlist
        j=j+1
    print("Chain not stabilized within MAXITER iterations")
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
    return False


