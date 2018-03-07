function main_fmm

clear
rand( 'seed',0)
randn('seed',0)

%%%%%%%%%%%%%%% Set the geometry.

  nmax = 50;
  ntot = 900;
  tt = 2*pi*rand(1,ntot);
  rr = 1 + 0.025*randn(1,ntot);
  xx = [rr.*cos(tt);...
        rr.*sin(tt)];
  x1min = min(xx(1,:));
  x1max = max(xx(1,:));
  x2min = min(xx(2,:));
  x2max = max(xx(2,:));
  len   = (1 + 1e-10)*max(x1max - x1min,x2max - x2min);
  box_geom_root = [len,0.5*(x1min+x1max),0.5*(x2min+x2max)];

%%%%%%%%%%%%%%%% Set various parameters
flag_precomp = 0;    % Determines whether to precompute the operators.
                     % (This applies _only_ to operators that depend on xx.)
flag_mode    = '11'; % This parameter specifies the type of the 
                     % charges and the sources.
                     %   flag_mode = '00'   monopoles given, potentials sought
                     %   flag_mode = '01'   monopoles given, fields sought
                     %   flag_mode = '10'   dipoles given,   potentials sought
                     %   flag_mode = '11'   dipoles given,   fields sought
p_fmm        = 30;   % The order of the multipole expansions used.
params       = ntot; % The total number of points.
fprintf(1,'======= FMM test code ==============\n')
fprintf(1,'ntot      = %d\n',ntot)
fprintf(1,'p_fmm     = %d\n',p_fmm)
fprintf(1,'nbox_max  = %d\n',nmax)
fprintf(1,'flag_mode = %s\n',flag_mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform the actual compression %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NODES,T_OPS,indvec] = ...
     LOCAL_FMM01_init(xx,nmax,p_fmm,flag_precomp,flag_mode,params);
nlevels = NODES{2,size(NODES,2)};                               

fprintf(1,'t_init    = %0.2e (sec)\n',t_init)
mem1 = whos('NODES');
mem2 = whos('T_OPS');
mem  = mem1.bytes + mem2.bytes;
fprintf(1,'memory    = %0.2e (MB)\n',mem/(2^20))
fprintf(1,'nlevels   = %d\n',nlevels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set up and solve a test problem %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_mode(2) == '0') % Monopole charges.
  qq = randn(ntot,3);
elseif (flag_mode(2) == '1') % Dipole charges.
  qq = randn(ntot,3) + 1i*randn(ntot,3);
end

uu = LOCAL_FMM01_apply(xx,NODES,T_OPS,indvec,qq,flag_precomp,flag_mode,params);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NODES,T_OPS,indvec] = LOCAL_FMM01_init(xx_orig,nmax,p_fmm,flag_precomp,flag_mode,params)

%%% Construct a root box.
x1max = max(max(xx_orig(1,:)));
x1min = min(min(xx_orig(1,:)));
x2max = max(max(xx_orig(2,:)));
x2min = min(min(xx_orig(2,:)));
len   = (1 + 1e-10)*max(x1max - x1min, x2max - x2min);
box_geom_root    = zeros(1,3);
box_geom_root(1) = len;
box_geom_root(2) = 0.5*(x1min + x1max);
box_geom_root(3) = 0.5*(x2min + x2max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Set up a tree structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NODES,nlevels,indvec] = LOCAL_get_tree(xx_orig,box_geom_root,nmax);
nboxes                 = size(NODES,2);
xx                     = xx_orig(:,indvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Create the interaction lists.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a list "leaf_list" of all childless boxes.
leaf_list  = zeros(1,nboxes);
ileaf      = 0;
for ibox = 1:nboxes;
  if (length(NODES{04,ibox}) == 0)
    ileaf = ileaf+1;
    leaf_list(ileaf) = ibox;
  end
end
nleaf = ileaf;
leaf_list = leaf_list(1:nleaf);

% Initialize all lists.
for ibox = 1:nboxes
  for ilist = 10:2:18
    NODES{ilist,ibox} = [];
  end
end

% Create lists 2 and 5
for ibox = 2:nboxes
  lside          = NODES{1,ibox}(1);
  x1c            = NODES{1,ibox}(2);
  x2c            = NODES{1,ibox}(3);
  ifath          = NODES{3,ibox};
  fath_coll_list = [NODES{18,ifath},ifath];
  for ifathcoll = fath_coll_list
    for jbox = NODES{04,ifathcoll}
      if ~(ibox == jbox)
        y1c  = NODES{1,jbox}(2);
        y2c  = NODES{1,jbox}(3);
        dist = sqrt((y1c-x1c)*(y1c-x1c) + (y2c-x2c)*(y2c-x2c));
        if (dist < (1.9*lside))
          % jbox is a colleague of ibox
          NODES{18,ibox} = [NODES{18,ibox},jbox];
        else
          NODES{12,ibox} = [NODES{12,ibox},jbox];
        end
      end
    end
  end
end

% Create lists 1 and 3
for ibox = 2:nboxes
  ibox_geom = NODES{1,ibox};
  list1   = [];
  list3   = [];
  if (length(NODES{04,ibox}) == 0)
    for jbox = NODES{18,ibox}
      [list1_fromjbox,list3_fromjbox] = LOCAL_get_list1_and_list3(NODES,ibox_geom,jbox,xx);
      list1 = [list1,list1_fromjbox];
      list3 = [list3,list3_fromjbox];
    end
  end
  NODES{10,ibox} = list1;
  NODES{14,ibox} = list3;
end

% At this point, NODES{10,ibox} lists only touching boxes at the same
% or finer levels as ibox. We next add the missing boxes.
for ibox = 2:nboxes
  for jbox = NODES{10,ibox}
    lside_ibox = NODES{1,ibox}(1);
    lside_jbox = NODES{1,jbox}(1);
    if ((lside_ibox / lside_jbox) > 1.0001)
      NODES{10,jbox} = [NODES{10,jbox},ibox];
    end
  end
end

% Create list 4.
% It is the dual of list 3.
for ibox = 2:nboxes
  for jbox = NODES{14,ibox}
    NODES{16,jbox} = [NODES{16,jbox},ibox];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [list1,list3] = LOCAL_get_list1_and_list3(NODES,ibox_geom,jbox,xx)

lside_i = ibox_geom(1);
lside_j = NODES{1,jbox}(1);

lsep    = max(abs(NODES{1,jbox}(2)-ibox_geom(2)),abs(NODES{1,jbox}(3)-ibox_geom(3)));
          
is_touching  = (lsep < 0.50001*(lside_i + lside_j));
has_children = (length(NODES{04,jbox})>0);

if ~is_touching
  list1 = [];
  list3 = [jbox];
else
  if ~has_children
    list1 = [jbox];
    list3 = [];
  else
    list1 = [];
    list3 = [];
    for kbox = NODES{04,jbox};
      [list1_fromk,list3_fromk] = LOCAL_get_list1_and_list3(NODES,ibox_geom,kbox,xx);
      list1 = [list1,list1_fromk];
      list3 = [list3,list3_fromk];
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NODES,nlevels,indvec] = LOCAL_get_tree(xx,box_geom,nmax)

ntot = size(xx,2);
len  = box_geom(1);
x1c  = box_geom(2);
x2c  = box_geom(3);


% BOXES is a temporary work array for setting up the tree-structure
BOXES       = zeros(18,100);
lenBOXES    = 100;
BOXES( 1,1) = NaN;
BOXES( 2,1) = 0;
BOXES( 3,1) = NaN;
BOXES( 6,1) = 1;
BOXES( 7,1) = ntot;
BOXES(10,1) = NaN;
BOXES(11,1) = x1c - 0.5*len;
BOXES(12,1) = x1c + 0.5*len;
BOXES(13,1) = x2c - 0.5*len;
BOXES(14,1) = x2c + 0.5*len;
BOXES(15,1) = -1;
BOXES(16,1) = -1;
BOXES(17,1) = -1;
BOXES(18,1) = -1;

indvec = 1:ntot;

% Create the tree structure by quartering any
% box that holds more than nmax nodes.
%
% We create the boxes one level at a time.
% The following loop is over LEVELS.
ibox_last  = 0;
ibox_new   = 1;
ilevel     = 0;
while (ibox_new > ibox_last)

  ibox_first = ibox_last+1;
  ibox_last  = ibox_new;
  ilevel     = ilevel+1;
  
  % Loop over all boxes on the level that was last created.
  % All newly created boxes are temporarily stored in the array TMPBOXES.
  for ibox = ibox_first:ibox_last
    
    % If ibox holds more than nmax nodes, it will be partitioned.
    if (BOXES(07,ibox) > nmax) 

      qfirst = BOXES(06,ibox);
      x1min  = BOXES(11,ibox);
      x1max  = BOXES(12,ibox);
      x2min  = BOXES(13,ibox);
      x2max  = BOXES(14,ibox);
      x1half = (x1min + x1max)/2;
      x2half = (x2min + x2max)/2;
      J      = qfirst - 1 + (1:BOXES(07,ibox));
      indloc = indvec(J);
      J_sw   = find( (xx(1,indloc) <= x1half) .* (xx(2,indloc) <= x2half) );
      J_nw   = find( (xx(1,indloc) <= x1half) .* (xx(2,indloc) >  x2half) );
      J_se   = find( (xx(1,indloc) >  x1half) .* (xx(2,indloc) <= x2half) );
      J_ne   = find( (xx(1,indloc) >  x1half) .* (xx(2,indloc) >  x2half) );
           
      indvec(J) = [indloc(J_sw),indloc(J_nw),indloc(J_se),indloc(J_ne)];
      
      npart_added = 0;
      % If there is not enough space to save the 4 children in
      % the array BOXES, then double the size of BOXES.
      if ((ibox_new + 4) > lenBOXES)
        BOXES  = [BOXES,zeros(size(BOXES,1),4+size(BOXES,2))];
        lenBOXES = size(BOXES,2);
      end
      if ~isempty(J_sw)
        ibox_new = ibox_new + 1;
        n_sw = length(J_sw);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_sw,NaN,NaN,NaN,x1min,x1half,x2min,x2half,...
                             -1,-1,-1,-1]';
        BOXES(15,ibox) = ibox_new;
        npart_added = npart_added + n_sw;
      end
      if ~isempty(J_nw)
        ibox_new = ibox_new + 1;
        n_nw = length(J_nw);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_nw,NaN,NaN,NaN,x1min,x1half,x2half,x2max,...
                             -1,-1,-1,-1]';
        BOXES(16,ibox) = ibox_new;
        npart_added = npart_added + n_nw;
      end
      if ~isempty(J_se)
        ibox_new = ibox_new + 1;
        n_se = length(J_se);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_se,NaN,NaN,NaN,x1half,x1max,x2min,x2half,...
                             -1,-1,-1,-1]';

        BOXES(17,ibox) = ibox_new;
        npart_added = npart_added + n_se;
      end
      if ~isempty(J_ne)
        ibox_new = ibox_new + 1;
        n_ne = length(J_ne);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_ne,NaN,NaN,NaN,x1half,x1max,x2half,x2max,...
                             -1,-1,-1,-1]';
        BOXES(18,ibox) = ibox_new;
        npart_added = npart_added + n_ne;
      end
    
    end
    
  end
  
end

nlevels = ilevel - 1;

% Let nboxes denote the number of boxes created,
% create the object NODES which will hold all
% relevant information, and transfer the 
% information in BOXES to NODES.
% We also delete the object BOXES.
nboxes  = ibox_new;
NODES = cell(45,nboxes);
for ibox = 1:nboxes
  NODES{01,ibox} = [     BOXES(12,ibox)-BOXES(11,ibox),...
                    0.5*(BOXES(11,ibox)+BOXES(12,ibox)),...
                    0.5*(BOXES(13,ibox)+BOXES(14,ibox))];
  NODES{02,ibox} = BOXES(02,ibox);
  NODES{03,ibox} = BOXES(03,ibox);
  NODES{04,ibox} = [];
  for j = 15:18
    if (BOXES(j,ibox) > 0)
      NODES{04,ibox} = [NODES{04,ibox},BOXES(j,ibox)];
    end
  end
  NODES{06,ibox} = BOXES(06,ibox);
  NODES{07,ibox} = BOXES(07,ibox);
end


return