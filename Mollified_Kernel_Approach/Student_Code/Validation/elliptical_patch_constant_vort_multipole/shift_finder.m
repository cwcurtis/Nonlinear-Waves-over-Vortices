function [dshift,ushift,lshift,rshift,ctgry] = shift_finder(col,row,rcnt,ccnt)
    if row > 1 && row < rcnt && col > 1 && col < ccnt        
        dshift = 1;
        ushift = 1;
        lshift = 1;
        rshift = 1;
        ctgry = 'interior';               
    elseif row == 1        
        dshift = 0;
        ushift = 1;            
        if col == 1 
            lshift = 0;
            rshift = 1;
            ctgry = 'tlcorner';            
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
            ctgry = 'trcorner';            
        else
            lshift = 1;
            rshift = 1;
            ctgry = 'thedge';            
        end        
    elseif row == rcnt    
        dshift = 1;
        ushift = 0;
        if col == 1
            lshift = 0;
            rshift = 1;
            ctgry = 'blcorner';            
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
            ctgry = 'brcorner';
        else
            lshift = 1;
            rshift = 1;
            ctgry = 'bhedge';
        end        
    else
        dshift = 1;
        ushift = 1;
        if col == 1
            lshift = 0;
            rshift = 1;
            ctgry = 'lvedge';        
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
            ctgry = 'rvedge';        
        end
    end