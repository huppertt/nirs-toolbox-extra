% meas = defiber(data, x1, x2, thresh, debug);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[meas] = defiber(data, x1, x2, thresh, debug)

if (~exist('debug','var') | isempty(debug))
   debug = 0;
end

h = figure;
set(h, 'doublebuff','on');

% I need a way to keep track of which pixels I've visited.  I can either
% use a global variable or pass something back with every function call.
% Of the two, the global variable is simplier.

global MASK;

MASK = zeros(size(data));

imagesc(MASK(620:720,520:620),[0 1]); drawnow;

% Check that we're within bounds

[n1,n2] = size(data);

if (debug)
   disp(sprintf('data is %d by %d', n1, n2));
end

if ((x1 < 1) | (x1 > n1))
   error('x1 out of range');
end

if ((x2 < 1) | (x2 > n2))
   error('x2 out of range');
end

% If the starting point is below threshold, then it's trivial

if (data(x1,x2) < thresh)
   if (debug)
      warning(sprintf('data(%d,%d) below threshold', x1, x2));
   end
   
   meas = 0;
else
   n0 = 1;
   y0 = data(x1,x2);
   
   [n1,y1] = bndfill(data, x1+1, x2, thresh, debug);
   [n2,y2] = bndfill(data, x1-1, x2, thresh, debug);
   [n3,y3] = bndfill(data, x1, x2+1, thresh, debug);
   [n4,y4] = bndfill(data, x1, x2-1, thresh, debug);
   
   meas = (y0 + y1 + y2 + y3 + y4) / (n0 + n1 + n2 + n3 + n4);
end

clear global MASK;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The real flood-fill algorithm

function[n,y] = bndfill(data, x1, x2, thresh, debug)

global MASK;

if (size(MASK) ~= size(data))
   warning('Size mismatch between mask and data');
   keyboard
end

if (debug)
   disp(sprintf('filling from data(%d,%d)',x1,x2));
end

% Check that we're within bounds

[m1,m2] = size(data);

n = 0;
y = 0;

if     ((x1 < 1) | (x1 > m1))
   if (debug)
      disp('-- x1 out of range');
   end  
elseif ((x2 < 1) | (x2 > m2))
   if (debug)
      disp('-- x2 out of range');
   end
elseif (MASK(x1,x2) ~= 0)
   if (debug)
      % Already visited
      disp(sprintf('-- data(%d,%d) masked out', x1, x2));
   end
elseif (data(x1,x2) < thresh)
   if (debug)
      % If the starting point is below threshold, then stop
      disp(sprintf('-- data(%d,%d) below threshold', x1, x2));
   end
else
   % Add in this point as well
   
   MASK(x1,x2) = MASK(x1,x2) + 1;
   imagesc(MASK(620:720,520:620),[0 1]); drawnow;
   
   [n1,y1] = bndfill(data, x1+1, x2, thresh, debug);
   [n2,y2] = bndfill(data, x1-1, x2, thresh, debug);
   [n3,y3] = bndfill(data, x1, x2+1, thresh, debug);
   [n4,y4] = bndfill(data, x1, x2-1, thresh, debug);
   
   n = n1 + n2 + n3 + n4 + 1;
   y = y1 + y2 + y3 + y4 + data(x1,x2);
   
   clear n1 n2 n3 n4 y1 y2 y3 y4;
end

disp(sprintf('Returns [%d,%f]', n, y));

return;
