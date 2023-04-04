import numpy

class Order:
    """
    Class implementing real quadratic orders.
    """
    def __init__(self, d: int):
        """
        Initializes a real qudratic order.
        If d = f^2*d_0 with d_0 squarefree, this initializes the unique order of conductor f in Q(sqrt(d)).
        Note that s = 1 if d_0 != 1 (mod 4) and = 2 otherwise, so that 
            omega_0 = (s - 1 + sqrt(d_0))/s 
        generates the ring of integers.
        """
        if(not is_square(d)):
            f = 1 
            d_0 = 1 
            for p in factor(d):
                e = p[1]
                while e > 1:
                    f *= p[0]
                    e -= 2
                if e == 1: d_0 *= p[0]
            if d_0%4 == 1: s = 2
            else: s = 1
            self.d = d
            self.d_0 = d_0
            self.s = s
            self.conductor = f
            if d_0%4 == 1: self.discriminant = f^2*d_0
            else: self.discriminant = 4*f^2*d_0
        else: self.d = self.d_0 = self.s = self.conductor = self.discriminant = 1
            
    def __str__(self):
        """
        Returns a string representation of the real quadratic order self.
        Output: (d = self.d, f = self.f, Delta = self.discriminant).
        Here d is the value used when initializing self, f is the conductor (d = f^2*d_0) and Delta the discriminant.
        Note that any of those 3 values uniquely determines the order.
        """
        return "(d = " + str(self.d) + ", f = " + str(self.conductor) + ", Delta = " + str(self.discriminant) + ")"
            
    def identity(self):
        """
        Returns the Z-module self = Z + Z*omega, implemented as the class Ideal.
        Here omega = f*omega_0 is uniquely determined.
        """
        return Ideal(self, self.s, self.conductor*(self.s-1))

    def regulator_cf(self, table = False, exact = 1, fu = False): 
        """
        Computes and returns the regulator of a real quadratic order, using the continued fraction algorithm.
        Settings: 
            - table: Controls whether a table of the principal cycle is also printed.
            - exact: Implements different levels of exactness for regulator computations:
                0: Fast (but potentially bad) approximations everywhere.  
                All outputs are exact. (Not recommended when table == True)
                1: The return value is a very good numerical approximation of the exact values,
                    all other distances are approximated (faster, but potentially worse).
                2: All outputs are very good numerical approximations of the exact values
                3: All outputs are exact. (Not recommended when table == True)
            - fu: Controls whether the fundamental unit (x + y*sqrt(d))/2 is also returned.
                    If True is passed, it is returned as a triple (x, y, o), where o is the norm.
                    This is the minimal nontrivial solution to x^2 - d*y^2 = +-4
                 
        """
        sqrtd_ = int(sqrt(self.d))                 
        if(exact == 0 or exact == 1): sqrtd = numpy.sqrt(self.d)        
        else: sqrtd = sqrt(self.d)
        a = self.identity()
        Q = Q_0 = Q_ = a.Q
        P = a.P
        q = (P + sqrtd_)//Q
        B_0, B_1, G_0, G_1 = 0, 1, Q, Q*q - P
        n = 1
        if table: print(n, Ideal(self, Q, P), 0)
        P = q*Q-P
        Q = (self.d - P^2)//Q
        q = (P+sqrtd_)//Q
        swap = G_1
        G_0, G_1 = swap, q*G_1 + G_0
        swap = B_1
        B_0, B_1 = swap, q*B_1+B_0
        if(exact == 0 or exact == 1): d = numpy.log((P + sqrtd)/Q_)
        else: 
            psi = (P + sqrtd)/Q_
            theta = psi
        while(Q != Q_0):
            n += 1
            if table:
                if(exact == 0 or exact == 1): print(n, Ideal(self, Q, P), d) 
                elif(exact == 2): print(n, Ideal(self, Q, P), numerical_approx(log(theta)))
                else: print(n, Ideal(self, Q, P), log(theta))
            Q_ = Q
            P = q*Q - P
            Q = (self.d - P^2)//Q
            q = (P + sqrtd_)//Q
            swap = G_1
            G_0, G_1 = swap, q*G_1 + G_0
            swap = B_1
            B_0, B_1 = swap, q*B_1 + B_0
            if(exact == 0 or exact == 1): d += numpy.log((P + sqrtd)/Q_)
            else: 
                psi = (P + sqrtd)/Q_
                theta *= psi
        n += 1
        if table:
            if(exact == 0 or exact == 1): print(n, Ideal(self, Q, P), d) 
            elif(exact == 2): print(n, Ideal(self, Q, P), numerical_approx(log(theta)))
            else: print(n, Ideal(self, Q, P), log(theta))
        if(exact == 0): ret = d 
        elif(exact == 3): ret = log((G_0 + sqrtd*B_0)/Q_0)
        else: ret = numerical_approx(log((G_0 + sqrtd*B_0)/Q_0))
        if fu: return (ret, (2*G_0//Q_0, 2*B_0//Q_0, (-1)^(n-1)))
        else: return ret
    
    def regulator_shanks(self, table = False, exact = 0):
        """
        Computes and returns the regulator of a real qudratic order, using Shanks' algorithm based on infrastructure.
        Settings:
                    - table: Controls whether a table of baby and giant steps performed is also printed.
                                This only happens if the algorithm does not return early.
                    - exact: Implements different levels of exactness:
                        -1: Using (f, p) representations for the giant steps (untested).
                        0: Distances are stored as numerical approximations. (No control over error propagation)
                        1: Keeps track of exact relative generators (and thus distances);
                            returns very good numerical approximations.
                        2: Keeps track of exact relative generators (and thus distances);
                            returns exact values. (Not recommended when table = True)
        """
        sqrtd_ = int(sqrt(self.d)) 
        if(exact == 0): 
            sqrtd = numpy.sqrt(self.d)
            sqrt4d = numpy.sqrt(sqrtd)
        else: 
            sqrtd = sqrt(self.d)
            sqrt4d = sqrt(sqrtd)
        a = self.identity()
        L = {}
        Q = Q_0 = Q_ = a.Q
        P = P_ = a.P
        q = (P + sqrtd_)//Q
        B_0, B_1, G_0, G_1 = 0, 1, Q, Q*q - P
        L[str(Ideal(self, Q, P))] = 0
        P = q*Q-P
        Q = (self.d - P^2)//Q
        q = (P+sqrtd_)//Q
        swap = G_1
        G_0, G_1 = swap, q*G_1 + G_0
        swap = B_1
        B_0, B_1 = swap, q*B_1+B_0
        psi = (P + sqrtd)/Q_
        if(exact == 0): lt = numpy.log(((P + sqrtd)/Q_))
        else:
            theta = psi
            lt = log(theta)
        if(Q == Q_0):
            if(exact == 0): return lt
            elif(exact == 1): return numerical_approx(lt)
            else: return log(theta)
        while(lt <= sqrt4d):
            if(exact == 1): L[str(Ideal(self, Q, P))] = numerical_approx(lt)
            else: L[str(Ideal(self, Q, P))] = lt
            Q_, P_ = Q, P
            P = q*Q - P
            Q = (self.d - P^2)//Q
            q = (P + sqrtd_)//Q
            swap = G_1
            G_0, G_1 = swap, q*G_1 + G_0
            swap = B_1
            B_0, B_1 = swap, q*B_1 + B_0
            psi = (P + sqrtd)/Q_
            if(exact == 0): lt += numpy.log((P + sqrtd)/Q_)
            else:
                theta *= psi
                lt = log(theta)
            if(Q == Q_0): 
                if(exact == 1): return numerical_approx(lt)
                else: return lt
            if(P_ == P): 
                if(exact == 0): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0/Q_)
                elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0/Q_))
                else: return 2*(lt - log(psi)) + log(Q_0/Q_)
            if(Q_ == Q): 
                if(exact == 0 or exact == 1): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0*psi/Q_)
                elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0*psi/Q_))
                else: return 2*(lt - log(psi)) + log(Q_0*psi/Q_)
        if(exact == 1): L[str(Ideal(self, Q, P))] = numerical_approx(lt)
        else: L[str(Ideal(self, Q, P))] = lt
        
        b = Ideal(self, Q, P)
        lt_b = lt
        if(exact != 0): theta_b = theta
        Q_, P_ = Q, P
        P = q*Q - P
        Q = (self.d - P^2)//Q
        q = (P + sqrtd_)//Q
        swap = G_1
        G_0, G_1 = swap, q*G_1 + G_0
        swap = B_1
        B_0, B_1 = swap, q*B_1 + B_0
        psi = (P + sqrtd)/Q_
        if(exact == 0): lt += numpy.log((P + sqrtd)/Q_)
        else:
            theta *= psi
            lt = log(theta)
        if(Q == Q_0): 
            if(exact == 1): return numerical_approx(lt)
            else: return lt
        if(P_ == P): 
            if(exact == 0): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0/Q_)
            elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0/Q_))
            else: return 2*(lt - log(psi)) + log(Q_0/Q_)
        if(Q_ == Q): 
            if(exact == 0 or exact == 1): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0*psi/Q_)
            elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0*psi/Q_))
            else: return 2*(lt - log(psi)) + log(Q_0*psi/Q_)
        if(exact == 1): L[str(Ideal(self, Q, P))] = numerical_approx(lt)
        else: L[str(Ideal(self, Q, P))] = lt
        Q_, P_ = Q, P
        P = q*Q - P
        Q = (self.d - P^2)//Q
        q = (P + sqrtd_)//Q
        swap = G_1
        G_0, G_1 = swap, q*G_1 + G_0
        swap = B_1
        B_0, B_1 = swap, q*B_1 + B_0
        psi = (P + sqrtd)/Q_
        if(exact == 0): lt += numpy.log((P + sqrtd)/Q_)
        else:
            theta *= psi
            lt = log(theta)
        if(Q == Q_0): 
            if(exact == 1): return numerical_approx(lt)
            else: return lt
        if(P_ == P): 
            if(exact == 0): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0/Q_)
            elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0/Q_))
            else: return 2*(lt - log(psi)) + log(Q_0/Q_)
        if(Q_ == Q): 
            if(exact == 0 or exact == 1): return 2*(lt - numpy.log(psi)) + numpy.log(Q_0*psi/Q_)
            elif(exact == 1): return numerical_approx(2*(lt - log(psi)) + log(Q_0*psi/Q_))
            else: return 2*(lt - log(psi)) + log(Q_0*psi/Q_)
        if(exact == 1): L[str(Ideal(self, Q, P))] = numerical_approx(lt)
        else: L[str(Ideal(self, Q, P))] = lt
        
        if table:
            print("Baby steps:")
            for l in L: print(l, L[l])
                
        if(exact == -1):
            fp_b = fp_b_ = b.to_fp_Rep(theta_b, 10)
            if table:
                print("Giant steps (fp):")
                print(fp_b)
            fp_b = fp_b.numult(fp_b_)
            while(str(fp_b.b) not in L):
                if table: print(fp_b)
                fp_b = fp_b.numult(fp_b_)
            if table: print(fp_b)
            return numerical_approx((fp_b.k*log(2) - L[str(fp_b.b)]))
                   
        b_ = b
        if(exact == 0): lt = lt_b
        else: theta = theta_b
        if table:
            print("Giant steps: ")
            if(exact == 1): print(b, numerical_approx(lt_b))
            else: print(b, lt_b)
        (b, (A, B, C)) = b.nucomp(b_, True)
        if(exact == 0):
            lt = lt + lt_b - numpy.log(abs((A + B*sqrtd)/C))
        else:
            theta = theta * theta_b * abs(1/((A + B*sqrtd)/C))
            lt = log(theta)
        while(str(b) not in L):
            if table:
                if(exact == 1): print(b, numerical_approx(lt))
                else: print(b, lt)
            (b, (A, B, C)) = b.nucomp(b_, True)
            if(exact == 0):
                lt = lt + lt_b - numpy.log(abs((A + B*sqrtd)/C))
            else:
                theta = theta * theta_b * abs(1/((A + B*sqrtd)/C))
                lt = log(theta)
        if table:
            if(exact == 1): print(b, numerical_approx(lt))
            else: print(b, lt)
                
        if(exact == 1): return numerical_approx(lt - L[str(b)])
        else: return lt - L[str(b)]
                         
class Ideal:
    """
    Class implementing ideals of a real quadratic order.
    """
    def __init__(self, O: Order, Q: int, P: int, S: int = 1):
        """
        Initializes an ideal.
        An ideal a is represented with generators as a Z-module:
            a = S*((Q/self.O.s)*Z + ((P + sqrt(d))/O.s)*Z),
        where O is the underlying order.
        """
        self.Q = Q
        self.P = P
        self.S = S
        self.O = O
        
    def __repr__(self):
        """
        Returns a string representation of an ideal a:
            a = (S)(Q, P) if S != 1 or
            a = (Q, P) if S = 1.
        """
        if self.S == 1: return "(" + str(self.Q) + "," + str(self.P) + ")"
        else: return "(" + str(self.S) + ")(" + str(self.Q) + "," + str(self.P) + ")"
    
    def reduce(self):
        """
        Reduces an ideal by applying the reduction step rho() until the norm guarantees that it is reduced.
        Untested and will only work in some cases (see the relevant theorems).
        """
        while(self.norm() >= numpy.sqrt(self.O.d)/2): self = self.rho()
        return self
        
    def rho(self):
        """
        Performs one reduction step rho(), 
        which is equivalent to one step in the simple continued fraction expansion of
            (P + sqrt(d))/Q.
        Untested.
        """
        #if(self.S != 1): return
        q = (self.P + int(numpy.sqrt(self.O.d)))//self.Q
        P = q*self.Q-self.P
        Q = (self.O.d-P^2)//self.Q
        #if Q < 0: return Ideal(self.O, -Q, -P)
        return Ideal(self.O, Q, P)
    
    def norm(self):
        """
        Returns the norm of an ideal
        """
        return self.S^2*self.Q//self.O.s
    
    def to_fp_Rep(self, theta: float, p):
        """
        Returns an (1, p) representation (b, d, k) of the ideal a, where a = theta*b
        """
        k = int(log(theta, 2))
        q = (int(2^(2*p-k)*theta) + 1)//2^p + 1
        r = (int(2^(2*p-k)*theta) + 1) - 2^p*q
        return fp_Rep(1, p, self, q, k)
        
    def multiply(self, b: Ideal):
        """
        Multiplies two ideals.
        Untested.
        """
        if(self.O is not b.O): return
        if(self.S != 1 or b.S != 1): return
        a = self
        if b.Q > a.Q:
            a, b = b, a
        var('x')
        G = gcd(a.Q//self.O.s, b.Q//self.O.s)
        try:
            X = Integer(solve_mod((b.Q//self.O.s)*x == G, a.Q//self.O.s)[0][0])
        except IndexError:
            X = 0        
        (S, Y, Z) = xgcd((a.P + b.P)//O.s, G)
        R_ = (self.O.d - b.P^2)//b.Q
        U = Integer(mod(X*Z*(a.P - b.P) + Y*R_, a.Q//S))     
        R_0, R_1, C_0, C_1, i = a.Q//S, U, 0, -1, 0
        Q = a.Q*b.Q//(self.O.s*S^2)
        P = Integer(mod(b.P + U*b.Q//(self.O.s*S), Q))
        return Ideal(self.O, Q, P, S)
        
    def nucomp(self, b: Ideal, gen = False):
        """
        Computes a reduced representative c in the ideal class [self*b] by performing NUCOMP.
        Settings:
            - gen: Controls whether integers A, B, C are returned, such that
                    c = |(A + B*sqrt(d))/C|*(self*a),
                i.e. (A + B*sqrt(d))/C is a relative generator of self*a with respect to c.
        """
        if(self.O is not b.O): return
        if(self.S != 1 or b.S != 1): return
        a = self
        if b.Q > a.Q:
            a, b = b, a
        var('x')
        G = gcd(a.Q//self.O.s, b.Q//self.O.s)
        try:
            X = Integer(solve_mod((b.Q//self.O.s)*x == G, a.Q//self.O.s)[0][0])
        except IndexError:
            X = 0
        (S, Y, Z) = xgcd((a.P + b.P)//self.O.s, G)
        R_ = (self.O.d - b.P^2)//b.Q
        U = Integer(mod(X*Z*(a.P - b.P) + Y*R_, a.Q//S))       
        R_0, R_1, C_0, C_1, i, j = a.Q//S, U, 0, -1, 0, 1
        if(R_0 < int(numpy.sqrt(2*self.O.s)*numpy.power(self.O.d, 1/4))):
            #print("j = 0")
            j = 0
            Q = a.Q*b.Q//(self.O.s*S^2)
            P = Integer(mod(b.P + U*b.Q//(self.O.s*S), Q))
            if gen:
                B_0, B_1 = 1, 0
            #    return(Ideal(self.O, Q, P), (S*(Q*B_0 + P*B_1), -S*B_1, Q))
            #else: return Ideal(self.O, Q, P)
        else:
            while(R_1 > int(numpy.sqrt(2*self.O.s)*numpy.power(self.O.d, 1/4))):
                i += 1
                q = R_0//R_1
                swap = C_1
                C_0, C_1 = swap, C_0 - q*C_1
                swap = R_1
                R_0, R_1 = swap, R_0 - q*R_1
            M_1, M_2 = ((b.Q//(self.O.s*S))*R_1 + (a.P - b.P)*C_1)//(a.Q//S), ((a.P + b.P)*R_1 + self.O.s*S*R_*C_1)//(a.Q//S)
            Q = (-1)^(i+1) * (R_1*M_1 - C_1*M_2)
            P = ((b.Q//(self.O.s*S))*R_1 + Q*C_0)//C_1 - b.P
        Q_ = abs(Q)
        k = (int(numpy.sqrt(self.O.d)) - P)//Q_
        P_ = k*Q_ + P
        if gen and j != 0:
            B_0, B_1 = abs(C_0), abs(C_1)
        if(P_ + int(sqrt(self.O.d)) >= Q_):
            #print("j = 1")
            if gen: return(Ideal(self.O, Q_, P_), (S*(Q*B_0 + P*B_1), -S*B_1, Q))
            else: return Ideal(self.O, Q_, P_)
        else:
            q = (P + int(numpy.sqrt(self.O.d)))//Q_
            Q_old = Q
            P = q*Q_ - P
            Q = (self.O.d - P^2)//Q_
            Q_ = abs(Q)
            k = int((int(numpy.sqrt(self.O.d)) - P)/Q_)
            P_ = k*Q_ + P
            if gen:
                swap = B_1
                B_0, B_1 = swap, q*B_1 + B_0
        if(P_ + int(numpy.sqrt(self.O.d)) >= Q_):
            #print("j = 2")
            B_1 *= sign(Q)
            if gen:
                return(Ideal(self.O, Q_, P_), (S*(Q*B_0 + P*B_1), -S*B_1, Q))
            else: return Ideal(self.O, Q_, P_)
        else:
            #print("j = 3")
            Q_old_ = Q
            Q = Q_ = Q_old - Q_old_ + 2*P
            P = Q_old_ - P
            P_ = P - Q
            if gen:
                swap = B_1
                B_0, B_1 = swap, B_1 + B_0
                return(Ideal(self.O, Q_, P_), (S*(Q*B_0 + P*B_1), -S*B_1, Q))
            else:
                return Ideal(self.O, Q_, P_)
           
#class fp_Rep:
#    """
#    Class implementing (f, p) representations of Ideals.
#    """
#    def __init__(self, f: float, p: int, b: Ideal, d: int, k: int):
#        """
#        Initializes an (f, p) representation (b, d, k).
#        No tests whether this is valid are performed.
#        """
#        self.f = numpy.float(f)
#        self.p = p
#        self.b = b
#        self.d = d
#        self.k = k
#        
#    def __repr__(self):
#        """
#        Returns a string representation ((f, p), (str(b), d, k))
#        """
#        return "(" + str(self.f) + ", " + str(self.p) + "): (" + str(self.b) + ", " + str(self.d) + ", " + str(self.k) + ")"
#    
#    def verify(self, theta: float):
#        """
#        Given a relative generator theta, returns whether the representation is valid.
#        Untested.
#        """
#        if(self.d < 2^self.p or self.d > 2^(self.p+1)): return False
#        if(abs((2^(self.p-self.k)*theta)/self.d - 1) >= self.f/2^self.p): return False
#        return True
#    
#    def remove(self, T: int, C: int, s: int, p: int):
#        """
#        Used as part of NUMULT
#        """
#        e = round(2^(self.p+3-s)*abs(T/C))
#        t = 0
#        while(2^(t-1) > e/(8*self.d) or e/(8*self.d) >= 2^t): t += 1
#        d_ = int(2^(p+3+t)*(self.d/e) + 1)
#        k_ = self.k - t
#        return fp_Rep(self.f + 9/8, self.p, self.b, d_, k_)
#    
#    def numult(self, fp: fp_Rep, optional = False): # optional missing
#        """
#        Implements NUMULT on (f, p)  representations.
#        For details see the relevant paper(s).
#        Untested.
#        """
#        if self.p != fp.p: return
#        p = self.p
#        if (self.b.O is not fp.b.O): return
#        (b, (A, B, C)) = self.b.nucomp(fp.b, True)
#        if(self.d*fp.d <= 2^(2*p+1)): e, h = (self.d*fp.d)//2^p, self.k + fp.k
#        else: e, h = (self.d*fp.d)//2^(p+1), self.k + fp.k + 1
#        s = 0
#        while(2^s*b.Q <= 2^(p+4)*B): s += 1
#        T = 2^s*A+B*int(2^s*sqrt(b.O.d))
#        f = self.f + fp.f + 2^(-self.p)*self.f*fp.f + 17/8
#        res = fp_Rep(f, p, b, e, h)
#        res = res.remove(T, C, s, p)
#        return res
