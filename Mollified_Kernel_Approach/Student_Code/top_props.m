function [ccnt,rcnt,nblcks,jvals] = top_props(ctgry)

if strcmp(ctgry,'tlcorner')
    ccnt = 4;
    rcnt = 4;
    nblcks = 16;
    jvals = [1;2;5;6];   
elseif strcmp(ctgry,'trcorner')
    ccnt = 4;
    rcnt = 4;
    nblcks = 16;
    jvals = [3;4;7;8];   
elseif strcmp(ctgry,'blcorner')
    ccnt = 4;
    rcnt = 4;
    nblcks = 16;
    jvals = [9;10;13;14];   
elseif strcmp(ctgry,'brcorner')
    ccnt = 4;
    rcnt = 4;
    nblcks = 16;
    jvals = [11;12;15;16];   
elseif strcmp(ctgry,'thedge')
    ccnt = 6;
    rcnt = 4;
    nblcks = 24;
    jvals = [3;4;9;10];    
elseif strcmp(ctgry,'bhedge')
    ccnt = 6;
    rcnt = 4;
    nblcks = 24;
    jvals = [15;16;21;22];    
elseif strcmp(ctgry,'lvedge')
    ccnt = 4;
    rcnt = 6;
    nblcks = 24;
    jvals = [9;10;13;14];    
elseif strcmp(ctgry,'rvedge')
    ccnt = 4;
    rcnt = 6;
    nblcks = 24;
    jvals = [11;12;15;16];    
elseif strcmp(ctgry,'interior')
    ccnt = 6;
    rcnt = 6;
    nblcks = 36;
    jvals = [15;16;21;22];    
end