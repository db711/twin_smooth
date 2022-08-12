#include <LiDIA/prime_list.h>
#include <LiDIA/quadratic_order.h>
#include <set>

using namespace LiDIA;
using namespace std;

void standard_form(const quadratic_number_standard & q, bigint & a, bigint & b, bigint & d);
/*
** Given a quadratic_number_standard computes the standard from
** 		((a + b*sqrt(disc))/2)/d
** that is required by some of our algorithms.
*/

bigint pp_coeff_modulo_bad(const quadratic_number_power_product & pp, bigint & x, bigint & y, const bigint & n);
/*
** Computes the coefficients x, y of a quadratic_number_power_product modulo some divisor of n (which is returned).
** The coefficients are computed with respect to the standard_form.
** This is Algorithm 2.
** Requires a normalized power product representation.
*/

void pp_coeff_modulo_good(const quadratic_number_power_product & pp, bigint & x, bigint & y, const bigint & n);
/* 
** Computes the coefficients x, y of a quadratic_number_power_product modulo.
** The coefficients are computed with respect to the standard_form.
** This is Algorithm 3.
** Requires a normalized power product representation.
*/

int pp_norm_sign(const quadratic_number_power_product & pp);
/*
** Computes the sign (+1 or -1) of the quadratic number representend in power product representation.
** Requires a normalized power product representation.
*/

int solve_special(const quadratic_number_power_product & pp);
/* 
** Computes and returns which power k of the fundamental unit corresponds to the solution of the special Pell equation.
** Operates according to Table 3.1.
*/

int pp_p_adic_valuation_y(const quadratic_number_power_product & pp, int p);
/*
** Computes the p-adic valuation of the y coefficient of the quadratic number given in power product representation.
** The computation is with respect to the standard form.
** Requires a normalized power product representation.
*/

bool pp_B_smooth_part_y(quadratic_number_power_product & pp, bigint & z, const prime_list & pl);
/*
** Given a list of prime numbers below B, this computes the part of the y coefficient of a qudratic number,
** given in power product representation, that is only divisible by those primes.
** This computation is NOT with respect to the standard form and instead assumed pp to be of the form
** 		x + y*sqrt(d)
** where d is the squarefree part of the discriminant.
** Requires a normalized power product representation.
**/

void pp_normalize(quadratic_number_power_product & pp);
/*
** Unlike the LiDIA documentation claims, power product representations are not always fully intialized,
** that is some if some factor in the product is the identity, it is omitted.
** As our algorithms require those factors to be present, we add them.
** Has to be called prior to most other functions.
*/

void twin_smooth(const string & input, const string & output, int B);
/*
** Given a list of discriminants and corresponding regulators, computes the corresponding
** twin smooth integers, according to Algorithm 1 (respectively Algorithm 5).
*/

int main()
{
    //twin_smooth("./input.txt", "./output.txt", B);
    return 0;
}

void twin_smooth(const string & input, const string & output, int B)
{
    ifstream in;
    ofstream out;
    const char * c;
    char * c_;
    string line;
    bigint D, d;
    bigfloat R;
    quadratic_order O;
    quadratic_number_power_product pp;
    bigint x, y, r, z = 1;
    prime_list pl;
    int k;
    bool prim_div;
    set<int> indices;

    pl.set_upper_bound(B);
    in.open(input);
    out.open(output);
    if(!out.is_open())
    {
		cout << "Error opening output file!" << endl;
        in.close();
        out.close();
        return;
    }
    if(in.is_open())
    {
        getline(in, line); getline(in, line); //skip first two lines
        while(getline(in, line, ' '))
        {
            c = line.c_str();
            string_to_bigint(c, D);
            getline(in, line);
            c_ = line.data(); //for some reason string_to_bigfloat needs writeable char * ...
            string_to_bigfloat(c_, R);
            O = quadratic_order(D);
            pp.assign(quadratic_ideal(O), R, 1); //assigns a compact representation of the fundamental unit
            pp_normalize(pp);
            remainder(r, D, 4);
            if(r == 0) d = D/4;
            else d = D;
            k = solve_special(pp);
            if(k > 1)
            {
                pp.assign(quadratic_ideal(O), k*R, 1)
                pp_normalize(pp);
            } //pp is now the fundamental solution to the special Pell equation
            indices.clear();
            for(int i = 1; i <= pl.get_last_prime() + 1; i++)
            {
                if(indices.count(i)) continue;
                if(i > 1)
                {
                    pp.assign(quadratic_ideal(O), k*i*R, 1);
                    pp_normalize(pp);
                }
               prim_div = pp_B_smooth_part_y(pp, z, pl); 
               if(abs(log(bigfloat(z))-i*k*R+log(2)+log(sqrt(bigfloat(d))))<.5) // yi is B-smooth
               {
                    pp_coeff_modulo_good(pp, x, y, 4); //4 to account for standard representation
                    if(x==2)
                    {
                        x = (pp.evaluate().get_a() - 1)/2; //(x, x+1) are B-smooth
                        //cout << x << "\n";
                        out << x << "\n";
                    }
               }
               else if(i == 1) break;
               else
               {
                   for(int j = 2*i; j <= pl.get_last_prime() + 1; j+=j) indices.insert(j);
               }
               if(prim_div && i >= 2) //this could be better, but we don't have a factorization of the discriminant
               {
                   if(i >= 6) break;
                   else
                   {
                       i = 5;
                       continue;
                   }
               }
            }
            //out.flush();
        }
    }
    else cout << "Error opening input file!" << endl;
    in.close();
    out.close();
    return;
}

void pp_normalize(quadratic_number_power_product & pp)
{
    quadratic_number_power_product_basis B = pp.get_basis(), C;
    base_vector< bigint > E = pp.get_exponents();
    quadratic_number_standard one;

    one.assign_one();
    one.set_order(pp.get_order());
    C.set_basis(one);
    while(E[E.get_size() - 1] != 1)
    {
        B.concat(B, C);
        E.insert_at(E[E.get_size() - 1]/2, E.get_size());
    }
    pp.set_basis(B);
    pp.set_exponents(E);
    return;
}

bool pp_B_smooth_part_y(quadratic_number_power_product & pp, bigint & z, const prime_list & pl)
{
    bigint  pow, r;
    bool b = true;
    int k;

    z = 1;
    k = pp_p_adic_valuation_y(pp, 2);
    remainder(r, pp.get_order().discriminant(), 4);
    if(r == 1) k--;
    if(k <= 0) b = false;
    else
    {
        power(pow, 2, k);
        z *= pow;
    }

    for(int i = 1; i < pl.get_number_of_primes(); i++)
    {
       k = pp_p_adic_valuation_y(pp, pl[i]); 
       if(k == 0) b = false;
       else
       {
            power(pow, pl[i], k);
            z *= pow;
       }
    }
    return b;
}

int pp_p_adic_valuation_y(const quadratic_number_power_product & pp, int p)
{
    int k = 10, ko = 0, m;
    bigint pow;
    bigint x, y;

    power(pow, p, k);
    pp_coeff_modulo_good(pp, x, y, pow);
    while(y == 0) 
    {
        ko = k;
        k *= 2;
        power(pow, p, k);
        pp_coeff_modulo_good(pp, x, y, pow);
    }

    while(1)
    {
        m = int((k - ko)/2. + ko);
        power(pow, p, m);
        pp_coeff_modulo_good(pp, x, y, pow);
        if(y == 0) ko = m;
        else k = m;
        if(ko + 1 == k) break;
    }
    return ko;
}

int solve_special(const quadratic_number_power_product & pp)
{
    bigint D = pp.get_order().discriminant();
    bigint r, x, y, d;
    int sign, k;
    
    remainder(r, D, 4);
    sign = pp_norm_sign(pp);
    if(r != 1)
    {
        if(sign == 1) k = 1;
        else k = 2;
    }
    else
    {
        pp_coeff_modulo_good(pp, x, y, 2);
        if(y == 0)
        {
            if(sign == 1) k = 1;
            else k = 2;
        }
        else
        {
            remainder(r, D, 2);
            if(r == 0) k = 2;
            else
            {
                if(sign == 1) k = 3;
                else k = 6;
            }
        }
    }
    return k;
}

int pp_norm_sign(const quadratic_number_power_product & pp)
{
    quadratic_number_power_product_basis B = pp.get_basis();
    bigint a, b, d;
    standard_form(B[B.get_size() - 1], a, b, d);
    remainder(a, a, 3);
    remainder(b, b, 3);
    remainder(d, B.get_order().discriminant(), 3);
    d = a*a - d*b*b;
    remainder(d, d, 3);
    if(d < 0) d += 3;
    if(d == 1) return 1; 
    else return -1;
}

void pp_coeff_modulo_good(const quadratic_number_power_product & pp, bigint & x, bigint & y, const bigint & n)
{
    bigint n_, m, k = 0;
    n_ = pp_coeff_modulo_bad(pp, x, y, n);
    while(n_ < n)
    {
        k += 10;
        power(m, n, k);
        n_ = pp_coeff_modulo_bad(pp, x, y, m);
        if(n_ < n) k += 10;
    }
    if(n < n_) pp_coeff_modulo_bad(pp, x, y, (m*n)/n_);
    return;
}

bigint pp_coeff_modulo_bad(const quadratic_number_power_product & pp, bigint & x, bigint & y, const bigint & n)
{
    quadratic_number_power_product_basis B = pp.get_basis();
    base_vector< bigint > E = pp.get_exponents();
    bigint a, b, d, di, xi, yi, r = 1, ri, m, s;
    bigint xold, yold;
    bigint D = 1;

    for(int i = 0; i < B.get_size() - 1; i++)
    {
        standard_form(B[i], a, b, d);
        D *= d*d;
    }
    standard_form(B[B.get_size() - 1], a, b, d);
    D *= d;
    power(d, 2, 2*B.get_size() - 1);
    D *= d*n;

    standard_form(B[0], a, b, d);
    remainder(x, a, D);
    remainder(y, b, D);

    for(int i = 1; i < B.get_size(); i++)
    {
        remainder(xi, (x*x + pp.get_order().discriminant()*y*y), D);
        remainder(yi, 2*x*y, D);
        standard_form(B[i], a, b, d);
        s = xi;
        remainder(xi, (xi*a + pp.get_order().discriminant()*yi*b), D);
        remainder(yi, (s*b + yi*a), D);
        r = r*r;
        standard_form(B[i-1], a, b, di);
        d = gcd(4*di*di*r, gcd(xi, yi));
        r = (4*di*di*r)/d;
        D = D/gcd(D, d);
        remainder(x, xi/d, D);
        remainder(y, yi/d, D);
        ri = r/gcd(r, D);
        while(gcd(ri, D) > 1) ri = ri/(gcd(ri, D));
        if(ri != 1)
        {
            power_mod(m, ri, -1, D);
            r = r/ri;
            remainder(x, x*m, D);
            remainder(y, y*m, D);
        }
        //cout << "i: " << i << "; " << x << ", " << y << " remainder " << r << " mod " << D << endl; 
    }
    standard_form(B[B.get_size() - 1], a, b, di);
    D = D/gcd(D, 2*di);
    remainder(x, x/di, D);
    remainder(y, y/di, D);
    if(x < 0) x += D; //remainders might be negative
    if(y < 0) y += D; //we prefer them to be positive
    return D;
}

void standard_form(const quadratic_number_standard & q, bigint & a, bigint & b, bigint & d)
{
    if(remainder(q.get_d(), 2))
    {
        a = 2*q.get_a();
        b = 2*q.get_b();
        d = q.get_d();
    }
    else
    {
        a = q.get_a();
        b = q.get_b();
        d = q.get_d();
        d.divide_by_2();
    }
    return;
}
