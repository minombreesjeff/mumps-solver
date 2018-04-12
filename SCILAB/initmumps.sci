function id = initmumps()
//
// id = initmumps
// it returns a default Scilab MUMPS mlist (structure)
//

id = mlist(["nom";"SYM";"JOB";"ICNTL";"CNTL";"PERM_IN";"COLSCA";"ROWSCA";"RHS";"INFOG";"RINFOG";"VAR_SCHUR";"SCHUR";"INST";"SOL";"DEFICIENCY";"NULLSPACE";"SYM_PERM";"UNS_PERM";"TYPE"],0,-1,zeros(1,40)-9998,zeros(1,5)-9998,-9999,-9999,-9999,-9999,zeros(1,40)-9998,zeros(1,20)-9998,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,0);

endfunction

