function lmisys = dlmiterm(lmisys, termID, A, B, flag, derivative)
%DLMITERM Add a term to a given group of LMIs representing a DLMI
%  
%     dlmiterm(LMISYS,TERMID,A,B,FLAG,IS_DERIVATIVE) adds one term to some
%     DLMI in the LMI system currently specified (see SETLMIS). A term is
%     either an outer factor, a constant matrix, or a variable term  A*X*B
%     or A*X'*B  where X is a matrix variable that can have a time
%     dependency or not. The first input TERMID specifies which LMI and
%     block the term belongs to.
%  
%     Note: Because the OFF-DIAGONAL blocks (i,j) and (j,i) are transposed of
%     ****  one another, specify the term content of ONLY ONE of these two 
%           blocks.
%  
%     Input:
%      LMISYS      current updated LMI system given by GETLMIS
%
%      TERMID      4-entry vector specifying the term location and nature
%                  Which DLMI?
%  		             TERMID(1) = +n  ->  left-hand side of the n-th LMI
%  		             TERMID(1) = -n  ->  right-hand side of the n-th LMI
%                  
%                    Note: Every DLMI will generate multiple LMIs
%                    ****  corresponding to the size of the DLMI variable
%                          (see DLMIVAR)
%                Which block?
%  		             For outer factors, set  TERMID(2:3) = [0 0].
%  		             Otherwise, set  TERMID(2:3) = [i j]  if the term
%                      belongs to the (i,j) block of the DLMI
%                What type of term?
%  		             TERMID(4) is an integer -> will be understood as the
%                        size of the DLMI, i.e., the number of LMIs used by
%                        the method to describe the DLMI
%                    TERMID(4) is a cell -> the cell is treated as a DLMI
%                        variable and the size of the DLMI will be derived
%                        by the size of the DLMI variable (see dlmivar)
%
%      A           value of the outer factor, constant term, or left
%                  coefficient in variable terms A*X*B or A*X'*B
%                  where X is the variable identifier returned by LMIVAR
%
%      B           right coefficient in variable terms A*X*B or A*X'*B
%
%      FLAG        quick way of specifying the expression  A*X*B+B'*X'*A'
%                  in a DIAGONAL block.  Set  FLAG='s'  to specify it
%                  with only one lmiterm command
%
%      DERIVATIVE  boolean that specifies whether the variable is meant to
%                  be interpreted as itself or its derivative, i.e., when
%                  DERIVATIVE is set to true, take dX/dt instead of X. The
%                  value for derivative is approximated by piecewise linear
%                  numerical method
%  
%     See also  setlmis, lmivar, getlmis, lmiedit, newlmi.
% 
%     Documentation for lmiterm

    arguments
        lmisys; termID;
        A = 1; B = 1;
        flag (1, 1) char {} = 'f';
        derivative (1, 1) logical {} = false;
    end
    setlmis(lmisys);
    [X, size_X, termID, term_is_constant] = handle_termID(termID);
    tag = abs(termID{1});
    side = sign(termID{1});
    counter = 0;
    if flag == 's' && ~term_is_constant
        for i = 1:size_X-1
            if derivative
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B, flag);
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, -B, flag);
                counter = counter + 1;
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B, flag);
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, -B, flag);
                counter = counter + 1;
            else
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B, flag);
                counter = counter + 1;
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, B, flag);
                counter = counter + 1;
            end
        end
    elseif ~term_is_constant
        for i = 1:size_X-1
            if derivative
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B);
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, -B);
                counter = counter + 1;
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B);
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, -B);
                counter = counter + 1;
            else
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i}], A, B);
                counter = counter + 1;
                lmiterm([side*(tag+counter), termID{2}, termID{3}, X{i+1}], A, B);
                counter = counter + 1;
            end
        end
    else
        for i = 1:size_X-1
            lmiterm([side*(tag+counter), termID{2}, termID{3}, 0], A);
            counter = counter + 1;
            lmiterm([side*(tag+counter), termID{2}, termID{3}, 0], A);
            counter = counter + 1;
        end
    end
end

function [var, var_size, new_termID, term_is_constant] = handle_termID(termID)
    var = termID(4:end);
    if ~isa(var, 'cell')
        var_size = var;
        var = repmat({0}, 1, var_size);
        new_termID = num2cell(termID);
        term_is_constant = true;
    else
        new_termID = termID;
        term_is_constant = false;
        var_size = size(var, 2);
    end
end
