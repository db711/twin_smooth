import itertools
#import quadratic_orders

def stormer(B: int, b: int = -1):
    """
    Computes all twin B-Smooth integers, using the continued fraction method to solve Pell equations
    as defined in quadratic_orders.py.
    Starts computation with smaller discriminants.
    """
    S = prime_range(B+1)
    T = set()
    #T = {}
    for d in stormer_discriminants(B):
        if(b != -1 and d > b): break
        (x, y, o) = Order(d).regulator_cf(False, 0, True)[1]
        (x_0, y_0) = (x, y) = solve_special(x, y, o, d)
        I = set(range(1, S[-1] + 2))
        i = 1
        while(i in I):
            if(is_smooth(y, B)): 
                if(is_odd(x)): T.add((x-1)//2) 
                #if(is_odd(x)): T[(x-1)//2] = d 
            else:
                if(i == 1): break
                else: I.intersection(set(range(i, S[-1] + 2, i)))
            j = 1
            while(not (i+j) in I and j < S[-1] + 2): j += 1
            i += j
            if(i in I): (x, y) = advance_special(2*x, 2*y, x_0, y_0, j, d)
            # check if each primitive divisor already occured
    return T

def stormer_discriminants(B):
    """
    Given a smoothness bound B, let p_1, ..., p_n be all primes less than B.
    This function yields all elements of the set
        {\prod_{i=1}^n p_i^{e_i} | e_i \in {0, 1}}\{1}
    in ascending order.
    """
    S = prime_range(B+1)
    numbers = list(S)
    indices = [0] * len(S)
    R = [2]
    numbers[0] = 6
    indices[0] = 2
    yield(2)
    while(numbers[-1] != float("inf")):
        n = min(numbers)
        R.append(n)
        for i in range(len(S)):
            if(numbers[i] == n):
                j = indices[i]
                while(j < len(R) and R[j]%S[i] == 0): j += 1
                if(j < len(R)):
                    numbers[i] = R[j]*S[i]
                    indices[i] = j + 1
                    if(i == len(S) - 1):
                        R = R[indices[-1]:]
                        for k in range(len(S)):
                            indices[k] -= indices[-1]
                else: numbers[i] = float("inf")
        yield(n)
        
def random_stormer_discriminant(B, n):
    """
    Given a smoothness bound B, let p_1, ..., p_n be all primes less than B.
    This function returns a random element of
        {\prod_{i=1}^n p_i^{e_i} | e_i \in {0, 1}}\{1}
    having precisely n prime factors.
    The probability distribution is uniform.
    """
    S = prime_range(B+1)
    d = 1
    for k in range(n):
        i = randrange(0, len(S))
        while(d%S[i] == 0): i = randrange(0, len(S))
        d *= S[i]
    return(d)
                
def powerset(iterable):
    """
    Returns the powerset
    """
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1))

def solve_special(x_0: int, y_0: int, o_0: int, d: int): 
    """
    Given the coefficients of the fundamental unit corresponding to d,
    returns the solution to the special Pell equation.
    """
    (x, y, o) = (x_0, y_0, o_0) 
    if x%2==0:
        if y%2 == 0 and o == 1: k = 1
        else: k = 2
    else: 
        if o==1: k = 3
        else: k = 6
    (x,y) = advance_special(x, y, x_0, y_0, k-1, d)
    return(x//2, y//2)

def advance_special(x: int, y: int, x_0: int, y_0: int, j: int, d: int): 
    """
    Advances a solution to the special Pell (x, y) equation by k steps.
    Requires knowing the fundamental solution (x_9, y_0) to the special Pell equation.
    """
    if(j < 0): return
    if(j == 0): return(x, y)
    x_, y_ = 2, 0
    while(j > 1): # Compute (x_0 + y_0*sqrt(d))^j
        if(is_even(j)):
            swap = x_0
            x_0 = (x_0^2 + y_0^2*d)//2
            y_0 = y_0*swap
            j = j//2
        else:
            swap = x_
            x_ = (x_0*x_ + y_0*y_*d)//2
            y_ = (x_0*y_ + y_0*swap)//2
            swap = x_0
            x_0 = (x_0^2 + y_0^2*d)//2
            y_0 = y_0*swap
            j = (j-1)//2
    swap = x_0
    x_0 = (x_0*x_ + y_0*y_*d)//2
    y_0 = (swap*y_ + y_0*x_)//2
    swap = x
    x = (x*x_0 + y*y_0*d)//2
    y = (swap*y_0 + x_0*y)//2
    return(x, y)
        
def is_smooth(n: int, B: int):
    """
    Tests whether n is B-smooth.
    """
    for p in prime_range(B+1):
        n = n//(p^n.valuation(p))
    if n == 1: return True
    else: return False