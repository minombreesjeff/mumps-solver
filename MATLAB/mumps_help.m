%****************************************
%This help menu gives details about the use of dmumps, zmumps and initmumps
%****************************************
%
%--------------- Input Parameters ---------------
%
% - mat: sparse matrix which has to be provided as the second argument of dmumps if id.JOB is strictly larger than 0.
%
% - id.SYM: controls the matrix type (symmetric positive definite, symmetric indefinite or unsymmetric) and it has do be initialized by the user before the initialization phase of MUMPS (see id.JOB). Its value is set to 0 after the call of initmumps.
%
% - id.JOB: defines the action that will be realized by MUMPS: initialize, analyze and/or factorize and/or solve and release MUMPS internal C/Fortran data. It has to be set by the user before any call to MUMPS (except after a call to initmumps, which sets its value to -1).
%
% - id.ICNTL and id.CNTL: define control parameters that can be set after the initialization call (id.JOB = -1). See Section ``Control parameters'' of the MUMPS user's guide for more details. If the user does not modify an entry in id.ICNTL then MUMPS uses the default parameter. For example, if the user wants to use the AMD ordering, he/she should set id.ICNTL(7) = 0. Note that the following parameters are inhibited because they are automatically set within the interface: id.ICNTL(19) which controls the Schur complement option and id.ICNTL(20) which controls the format of the right-hand side.
%
% - id.PERM\_IN: corresponds to the given ordering option (see Section ``Input and output parameters'' of the MUMPS user's guide for more details). Note that this permutation is only accessed if the parameter id.ICNTL(7) is set to 1.
%
% - id.COLSCA and id.ROWSCA: are optional scaling arrays (see Section ``Input and output parameters'' of the MUMPS user's guide for more details)
%
% - id.RHS: defines the right-hand side. The parameter id.ICNTL(20) related to its format (sparse or dense) is automatically set within the interface. Note that id.RHS is not modified (as in MUMPS), the solution is returned in id.SOL.
%
% - id.VAR\_SCHUR: corresponds to the list of variables that appear in the Schur complement matrix (see Section ``Input and output parameters'' of the MUMPS user's guide for more details).
%
%--------------- Output Parameters ---------------
%
% - id.SCHUR: if id.VAR\_SCHUR is provided of size SIZE\_SCHUR, then id.SCHUR corresponds to a dense array of size (SIZE\_SCHUR,SIZE\_SCHUR) that holds the Schur complement matrix (see Section ``Input and output parameters'' of the MUMPS user's guide for more details). The user does not have to initialize it.
%
% - id.INFOG and id.RINFOG: information parameters (see Section ``Information parameters'' of the MUMPS user's guide ).
%
% - id.SYM\_PERM: corresponds to a symmetric permutation of the variables (see discussion regarding ICNTL(7) in Section ``Control parameters'' of the MUMPS user's guide ). This permutation is computed during the analysis and is followed by the numerical factorization except when numerical pivoting occurs.
%
% - id.UNS\_PERM: column permutation (if any) on exit from the analysis phase of MUMPS (see discussion regarding ICNTL(6) in Section ``Control parameters'' of the MUMPS user's guide ).
%
% - id.SOL: dense vector or matrix containing the solution after MUMPS solution phase.
%
%--------------- Internal Parameters ---------------
%
% - id.INST: (MUMPS reserved component) MUMPS internal parameter.
% 
% - id.TYPE: (MUMPS reserved component) defines the arithmetic (complex or double precision).
% 
