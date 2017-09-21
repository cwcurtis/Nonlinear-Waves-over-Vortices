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
            ctgry = 'corner';            
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
            ctgry = 'corner';            
        else
            lshift = 1;
            rshift = 1;
            ctgry = 'hedge';            
        end        
    elseif row == rcnt    
        dshift = 1;
        ushift = 0;
        if col == 1
            lshift = 0;
            rshift = 1;
            ctgry = 'corner';            
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
            ctgry = 'corner';
        else
            lshift = 1;
            rshift = 1;
            ctgry = 'hedge';
        end        
    else
        dshift = 1;
        ushift = 1;
        ctgry = 'vedge';
        if col == 1
            lshift = 0;
            rshift = 1;
        elseif col == ccnt
            lshift = 1;
            rshift = 0;
        end
    end