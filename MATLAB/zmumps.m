function [id]=zmumps(id,mat)
%
% [id]=zmumps(id,mat)
% id is a structure (see details in initmumps.m and MUMPS documentation)
% mat is optionnal if the job is -1 or -2
% mat is a square sparse matrice
% information are return in id fields
%
% Use help mumps_help for detailed information
%

errmsg = nargoutchk(1,1,nargout);
if(~isempty(errmsg))
     disp(errmsg);
     return;
end

prectype = 2;

if(id.JOB == -2)     
     if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
     end
     if(id.TYPE ~= prectype)
       disp('You are trying to call CMPLX/DBL version on a DBL/CMPLX instance');
       return;
     end
     zmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST);
     id = [];
     return;
end


if(id.JOB == -1)
     if(id.INST~=-9999)
         disp('Allready initialized instance');
         return;
     end
     [inform,rinform,sol,inst,schur,def,nullspace,sym_perm,uns_perm] = zmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST);
     id.INFOG = inform;
     id.RINFOG = rinform;
     id.SOL = sol;
     id.INST = inst;
     id.SCHUR = schur;
     id.DEFICIENCY = def;
     id.NULLSPACE = nullspace;
     id.SYM_PERM = sym_perm;
     id.UNS_PERM = uns_perm;
     id.TYPE = prectype;
     return;
end

if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
end

if(id.TYPE ~= prectype)
   disp('You are trying to call CMPLX/DBL version on a DBL/CMPLX instance');
   return;
end

[inform,rinform,sol,inst,schur,def,nullspace,sym_perm,uns_perm] = zmumpsmex(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,mat);
id.INFOG = inform;
id.RINFOG = rinform;
id.SOL = sol;
id.INST = inst;
if(id.JOB == 2 | id.JOB == 4 | id.JOB == 6)
  if(id.SYM == 0)
	id.SCHUR = schur';
  else
        id.SCHUR = triu(schur)+tril(schur',-1);
  end
end
id.DEFICIENCY = def;
id.NULLSPACE = nullspace;
id.SYM_PERM(sym_perm) = [1:size(mat,1)];
id.UNS_PERM = uns_perm;

