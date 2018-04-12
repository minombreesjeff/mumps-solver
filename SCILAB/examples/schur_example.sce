//A simple demo for the MUMPS interface, with the return of the schur complement
//to run it, You just have to execute  the instruction within Scilab
//		exec sparse_example.sce; 


//*********************** MATRIX INITIALISATION ***********************//
// This matrix has to be a SciSparse, otherwise it won't work.
	exec('ex2.sci');
	//voir pour les speyes
	mat=sparse(a);
	n = size(mat,1);
	
// Right Hand side setting
	RHS = ones(n,1);


	
//****************** Initialisation of the Scilab MUMPS structure ******************//
timer();
[id]=initmumps();

//Here Job=-1, the next call will only initialise the C and Fortran structure
[id]=dmumps(id);

id.RHS=RHS;
id.VAR_SCHUR = [n-4:n];
themax = max(max(abs(mat)));
mat = mat+sparse([1:n;1:n]',3*themax*ones(1,n));

//******************** CALL TO MUMPS FOR RESOLUTION ********************//
job=6;
id.JOB=job;

[id]=dmumps(id,mat);

// verification of the solution
solution=id.SOL;
norm1=norm(mat(1:n-5,1:n-5)*solution(1:n-5) - ones(n-5,1),'inf');
if norm1> 10^(-9) then
	write(%io(2),'WARNING: precision may not be OK');
else
	write(%io(2),'SCHUR SOLUTION: CHECK OK');
end

//****************** DESTRUCTION OF THE MUMPS INSTANCE ******************//
job=-2;
id.JOB=job;
[id]=dmumps(id);
t=timer()
