/*
** Computes the regulator of the maximal order in Q(sqrt(d)),
** where a list of values for d is given in an input file.
** Writes the discriminants and regulators into an outpute file
** and returns the time it took in seconds.
*/

regulators_file(input,output) = 
{
	local(m=fileopen(input,r),n=fileopen(output,w),d,start=getwalltime);
	\\default(parisize,4294967296); \\increase stacksize to 4GB
	fileread(m); \\skips first line
	while(d=eval(fileread(m)),
		if((d%4)!=1,d=4*d);
		if(d<10^10,
			filewrite1(n,d);filewrite1(n," ");filewrite(n,quadregulator(d)),
		\\else
			filewrite1(n,d);filewrite1(n," ");filewrite(n,quadclassunit(d)[4])
		)
	);
	fileflush(n);
	fileclose(n);
	fileclose(m);
	(getwalltime()-start)/1000.
}


/*
** Computes the regulators for all maximal orders in the real quadratic number field
** Q(sqrt(d)), where d is squarefree and from the set 
** {p_^(e_1)*...*p_n^(e_n) <= N | p_i <= B prime, e_i = 0, 1 for i = 1, ..., n}.
** Writes the discriminants and regulators to a file and returns the time it took in seconds.
** Setting N = -1 is treated as N = inf.
**/
stormer_regulators(B,N,file) =
{
	local(S=primes(primepi(B)),d=1,D,n=fileopen(file,w),start=getwalltime);
    \\default(parisize,4294967296); \\increase stacksize to 4GB
	forsubset(#S,s,
		D=vecextract(S,s);
		if(D==[],,
			d=1;
			foreach(D,p,d=d*p);
			if((d%4)!=1,d=4*d);
				if(N!=-1&&d>N,,
					if(d<10^10,
					filewrite1(n,d);filewrite1(n," ");filewrite(n,quadregulator(d)),
				\\else
					filewrite1(n,d);filewrite1(n," ");filewrite(n,quadclassunit(d)[4])
			)
			)
		)
	);
	fileflush(n);
	fileclose(n);
    (getwalltime()-start)/1000.
}

/*
** Down here are test functions to determine the best algorithm to use.
*/

\\ Return a random nonsqure between 10^n and 10^(n+1).
random_nonsquare(n) = 
{
	local(r=random(9^n)+10^n);
	while(issquare(r),r=random(9^n)+10^n);
	r
}

/*
** Computes 1000 regulators to randomly chosen discriminants 
** between 10^n and 10^(n+1) with the chosen method and returns time taken.
** t = 1, quadregulator(D)
** t = 2, quadclassunit(D)
** t = 3, bnfinit(x^2-d)
*/
test_regulator(n,t) = 
{
	local(start=getwalltime(),end,r=random_nonsquare(n));
	for(i=1,1000,
		if( t==1,if((r%4)!=1,r=4*r);quadregulator(r),
			t==2,if((r%4)!=1,r=4*r);quadclassunit(r)[4],
			t==3,bnfinit(x^2-r).reg,
			return(-1)
		);
		r=random_nonsquare(n)
	);
	end=getwalltime();
	(end-start)/1000.
}

/*
** Computes 1000 regulators to randomly chosen discriminants between
** 10^i and 10^(i+1) for each i=1,2,...,n with various methods:
** quadregulator(D), quadclassunit(D) and bnfinit(x^2-d).
** Prints timing data for each method.
**
** Tests show that quadregulator is best until 10^10,
** and quadclassunit beyond that.
*/
longtest_regulator(n) =
{
	for(i=1,n,
		print(i);
		for(j=2,3,
			if( j==1,print("quadregulator"),
				j==2,print("quadclassunit"),
				j==3,print("bnfinit")
			);
			print(test_regulator(i,j))
		);
		print()
	)
}
