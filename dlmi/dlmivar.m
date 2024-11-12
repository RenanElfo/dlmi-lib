function [X, n, sX] = dlmivar(type, structure, lmisys, size)
%DLMIVAR Declare dlmi variable as cell. Has similar structure to lmivar.
%
%   Create a new cell of given size of matrix-value variables in LMI system
%   which is interpreted as a dlmi variable when using dlmiterm.
%
%    X= DLMIVAR(TYPE,STRUCT,LMISYS,SIZE) adds a new dlmi variable X to
%    the LMI system specified by LMISYS. You can use the identifier X for
%    subsequent references to the variable X in calls to DLMITERM.
%     
%    [X,NDEC,XDEC] = DLMIVAR(TYPE,STRUCT,LMISYS,SIZE) also returns the
%    total number of decision variables associated with X, NDEC, and the
%    entry-wise dependence of X on these decision variables, XDEC.
%     
%        Input:
%         TYPE     Structure of X:
% 		      1 -> symmetric block diagonal
% 		      2 -> full rectangular
% 		      3 -> other
%         STRUCT   Additional data on the structure of X
% 	         TYPE=1: the i-th row of STRUCT describes the i-th diagonal block 
%                     of X
% 		             STRUCT(i,1) -> block size
% 		             STRUCT(i,2) -> block type, i.e., 
%                                    0 for scalar blocks t*I
% 			                        1 for full block
%                                    -1 for zero block
% 	         TYPE=2: STRUCT = [M,N]  if X is a MxN matrix
% 	         TYPE=3: STRUCT is a matrix of the same dimension as X
% 		         where STRUCT(i,j) is either
% 		           0  if X(i,j) = 0
% 		          +n  if X(i,j) = n-th decision variable
% 		          -n  if X(i,j) = (-1)* n-th decision var
%        LMISYS   The updated LMI system given by GETLMIS
%                    Must always use GETLMIS to get the updated LMI system.
%        SIZE     The size of the dlmi variable
%                    The given size of the dlmi variable must be coherent
%                    with the numerical method used. For instance,
%                    a piecewise linear numerical method must have dlmi
%                    variables of size 2^n+1 for a given n.
%        Output:
%         X        Optional: Identifier for the new matrix variable.
%                  Its value is k if k-1 matrix variables have already
%                  been declared.  This identifier is not affected by
%                  subsequent modifications of the LMI system
%         NDEC     total number of decision variables so far
%         XDEC     entry-wise dependence of X on the decision variables
%                  (Xdec = STRUCT for Type 3)
%     
%        Examples:
%         % Create a 2x2 full symmetric dlmi matrix variable of size 2^1+1
%         lmisys=getlmis;
%         size=2^1+1
%         X=dlmivar(1,[2 1],lmisys,size)
%     
%         % Create a 2x3 dlmi matrix variable Y of size 8
%         Y=dlmivar(2,[2 3],getlmis,8)
%     
%         % Create a 3x3 diagonal dlmi matrix variable Z of size 2
%         Z=dlmivar(1,[3 0],getlmis,2);
%
%   See also  lmivar, getlmis
   setlmis(lmisys);
   X = cell(1, size);
   sX = cell(1, size);
   for i = 1:size
       [X{i}, n, sX{i}] = lmivar(type, structure);
   end
end
