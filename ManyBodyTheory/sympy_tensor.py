# -*- coding: utf-8 -*-
"""
Use Sympy to get coupled cluster equations 
"""
"""
Calculates the Coupled-Cluster energy- and amplitude equations
See 'An Introduction to Coupled Cluster Theory' by
T. Daniel Crawford and Henry F. Schaefer III.
Other Resource : http://vergil.chemistry.gatech.edu/notes/sahan-cc-2010.pdf
"""

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, KroneckerDelta )
from sympy import (
    symbols, Rational, latex, Dummy
)

pretty_dummies_dict = {
    'above': 'cdefgh',
    'below': 'klmno',
    'general': 'pqrstu'
}


def get_CC_operators():
    """
    Returns a tuple (T1,T2) of unique operators.
    """
    i = symbols('i', below_fermi=True, cls=Dummy)
    a = symbols('a', above_fermi=True, cls=Dummy)
    t_ai = AntiSymmetricTensor('t', (a,), (i,))
    ai = NO(Fd(a)*F(i))
    i, j = symbols('i,j', below_fermi=True, cls=Dummy)
    a, b = symbols('a,b', above_fermi=True, cls=Dummy)
    t_abij = AntiSymmetricTensor('t', (a, b), (i, j))
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))
   
    T1 = t_ai*ai
    T2 = Rational(1, 4)*t_abij*abji
    return (T1, T2)

def get_T3_operator():
    """
    return T3 
    """
    i,j,k = symbols('i,j,k', below_fermi=True, cls=Dummy)
    a,b,c = symbols('a,b,c', above_fermi=True, cls=Dummy)
    t_abcijk = AntiSymmetricTensor('t', (a, b, c), (i, j, k))
    abckji = NO(Fd(a)*Fd(b)*Fd(c)*F(k)*F(j)*F(i))
    T3 = Rational(1,6)*t_abcijk*abckji
    return T3 

def lprint(x):
    # print in latex form 
    print(latex(x))
    return 

def main():
    print()
    print("Calculates the Coupled-Cluster energy- and amplitude equations")
    print("See 'An Introduction to Coupled Cluster Theory' by")
    print("T. Daniel Crawford and Henry F. Schaefer III")
    print("Reference to a Lecture Series: http://vergil.chemistry.gatech.edu/notes/sahan-cc-2010.pdf")
    print()

    #--------------------------
    # setup hamiltonian
    #
    #-------------------------
    p, q, r, s,t,u = symbols('p,q,r,s,t,u', cls=Dummy)
    f = AntiSymmetricTensor('f', (p,), (q,)) # f^p_q
    pr = NO(Fd(p)*F(q))    # a^+_p a_q 
    v = AntiSymmetricTensor('v', (p, q), (r, s)) # v^{pq}_{rs}
    pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r)) # a^+_p a^+_q a_s a_r 
    w = AntiSymmetricTensor('w', (p, q, r), (s, t, u)) # w^{pqr}_{stu}
    pqruts = NO(Fd(p)*Fd(q)*Fd(r)*F(u)*F(t)*F(s))

    H = f*pr + Rational(1, 4)*v*pqsr 
    
    Hp = f*pr + Rational(1, 4)*v*pqsr + Rational(1,6)*w*pqruts # this is for test 
    
    print("Using the hamiltonian:", latex(H))

    print("Calculating 4 nested commutators")
    C = Commutator

    T1, T2 = get_CC_operators() 
    T = T1 + T2
    print("commutator 1...")
    #--------------------
    # 1. C(H,T) just have [H,T] expression, not expanded yet
    # 2. wicks(C(H,T)) makes contractions using wicks theorem
    #               and replace contractions with  deltas 
    #               reults are list of Normal ordered operators. 
    #               no evaluation of delta is done. 
    #               thus quite lengthy 
    # 3. evaluate_deltas
    #          replace indices using deltas 
    #          f^{p}_{p} does not vanish though it is defined as 
    #                 anti-symmetric tensor. 
    #        one can use simplify_kronecker_deltas=True in wicks instead.
    # 4. substitute_dummies
    #          replace the dummy index to be easier, consistent form  
    #          resulting comm1 is list of normal ordered operators 
    #--------------------
    comm1 = wicks(C(H, T))
    comm1 = evaluate_deltas(comm1)
    comm1 = substitute_dummies(comm1)

    T1, T2 = get_CC_operators() 
    T = T1 + T2       # why define again ? seems to redundant 
    print("commutator 2...")
    comm2 = wicks(C(comm1, T))
    comm2 = evaluate_deltas(comm2)
    comm2 = substitute_dummies(comm2)

    T1, T2 = get_CC_operators()
    T = T1 + T2
    print("commutator 3...")
    comm3 = wicks(C(comm2, T))
    comm3 = evaluate_deltas(comm3)
    comm3 = substitute_dummies(comm3)

    T1, T2 = get_CC_operators()
    T = T1 + T2
    print("commutator 4...")
    comm4 = wicks(C(comm3, T))
    comm4 = evaluate_deltas(comm4)
    comm4 = substitute_dummies(comm4)

    print("construct Hausdorff expansion...")
    eq = H + comm1 + comm2/2 + comm3/6 + comm4/24
    eq = eq.expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, new_indices=True,
            pretty_indices=pretty_dummies_dict)
    print("*********************")
    print()

    print("extracting CC equations from full Hbar")
    i, j, k, l = symbols('i,j,k,l', below_fermi=True)
    a, b, c, d = symbols('a,b,c,d', above_fermi=True)
    print()
    print("CC Energy:")
    en=wicks(eq, simplify_dummies=True, keep_only_fully_contracted=True)
    print(latex(en))
    print()
    print("CC T1:")
    eqT1 = wicks(NO(Fd(i)*F(a))*eq, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
    eqT1 = substitute_dummies(eqT1)
    print(latex(eqT1))
    print()
    print("CC T2:")
    eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*eq, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    P = PermutationOperator
    eqT2 = simplify_index_permutations(eqT2, [P(a, b), P(i, j)])
    print(latex(eqT2))

    return en, eqT1, eqT2

#===Run===============
if __name__ == "__main__":
    
    en, eqT1, eqT2 = main()
    
    
    
    