function id=zmumps(id,mat)

//************************************************************************************************************** 
// [id] = zmumps(id,mat)
// id is a structure (see details in initmumps.m and MUMPS documentation)
// mat is an optional parameter if the job id.job = -1 or -2
// mat is a square sparse matrix
// informations are return in id fields
// 
// *************************************************************************************************************
prectype=1;

if id.JOB == -2 then   
     if id.INST==-9999 then
         disp('Uninitialized instance');
         return;
     end
     if id.TYPE ~= prectype then
	disp('You are trying to call CMPLX/DBL version on a DBL/CMPLX instance');
        return;
     end
     // call the C routine zmumpsc
 
     zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST);
     id = [];
     return;
end


if id.JOB == -1 then
	if id.INST~=-9999 then
		disp('Allready initialized instance');
	return;
	end
	// call the C routine zmumpsc
	[inform,rinform,sol,inst,schu,def,nullspace,sym_perm,uns_perm] = zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST);
        id.INFOG = inform;
	id.RINFOG = rinform;
	id.SOL = sol;
	id.INST = inst;
	id.SCHUR = schu;
	id.DEFICIENCY = def;
	id.NULLSPACE = nullspace;
	id.SYM_PERM = sym_perm;
	id.UNS_PERM = uns_perm;
        id.TYPE=prectype;

 	clear inform rinform sol inst schu def nullspace sym_perm uns_perm
	return;
	
end

if id.INST ==-9999 then
	disp('Uninitialized instance');
	return;
end 
// call the C routine zmumpsc

if id.TYPE ~= prectype then
	disp('You are trying to call CMPLX/DBL version on a DBL/CMPLX instance');
end

[inform,rinform,sol,inst,schu,def,nullspace,sym_perm,uns_perm] = zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,mat);
id.INFOG = inform;
id.RINFOG = rinform;
id.SOL = sol;
id.INST = inst;
if (id.JOB == 2|id.JOB==4|id.JOB==6) then
	if id.SYM == 0 then
                id.SCHUR=schu';
        else 
        	id.SCHUR=triu(schu)+tril(schu',-1);
        end
end 
id.DEFICIENCY = def;
id.NULLSPACE = nullspace;
id.SYM_PERM(sym_perm) = [1:size(mat,1)];
id.UNS_PERM = uns_perm;

clear inform rinform sol inst schu def nullspace sym_perm uns_perm

endfunction
