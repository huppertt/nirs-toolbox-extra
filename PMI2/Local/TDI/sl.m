% function[sl] = sl;

function[sl] = sl

% Source fibers 1, 3, 4, 16, 35, 48, 50 all seem to be broken.  Assign
%  holes to these fibers more or less at random since there's no good
%  way to find out where they actually go.

% These fibers are broken, they could go anywhere and I wouldn't be able to
% tell the difference.

holes = [ 1 3 16 50 ];

% Actual mapping

map = [ 55 41 35 45 24 47 39 ...
         5 25 42 26 28 14  0 ...
        57 51 30  0  7 12 56 ...
        31 23 32 19 17 59 38 ...
         4 13 48 20 60  0 49 ...
        33 29 46 62 22 21 40 ...
        11 52  9 61 44 53 63 ...
        36 58 37  2 54 18 15 ...
         8 10 43  6 27 34  0 ];

% Back-fill the holes in random order

ml = find(map == 0);

if (length(ml) ~= length(holes))
   warning(sprintf('hole/map mismatch [%d vs %d]', length(ml), length(holes)));
   keyboard
else
   % Backfill random points
   map(ml) = holes;
end

% Copy back

sl = map;

return;
