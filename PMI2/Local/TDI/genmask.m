% function[mask] = genmask(rawdata, threshold, doRight)
%
% rawdata   -> camera image used to generate mask [max(data,[],3) works best]
% threshold -> (rawdata > threshold) defines masked fibers
% doRight   -> right vs left detector positions [i.e., plug top/bottom]

function[mask] = genmask(rawdata, threshold, doRight)

% Pick which set the user should click on

if (~exist('doRight','var') | isempty(doRight))
   doRight = 1;
end

% Make sure I've got doubles to work with

rawdata = double(rawdata);

% Get user coordinates

h = figure;

imagesc(rawdata', [0 threshold]);

disp('Click on the 4 corners [UL, UR, LL, LR - IN THAT ORDER]')
disp('(not counting partial rows)');

r = ginput(4);

u = r(:,1);
v = r(:,2);

close(h);

% Mapping from recorded values to pixel coordinates

disp('Generating mask');

% Tabulate know spatial position of detector fibers

% "Right" detector fibers

XR = [ [  8  6  5  7  4  9  3 ] ...
       [  6  4  5  7  8 10  9 ] ...
       [  6 13 12 11 11 12 10 ] ...
       [  5  4  3  2  9  8 11 ]...
       [ 13 15 14 14 12 10  7 ] ...
       [  3  5  4  1  2  9 10 ] ...
       [ 13 11 12 12  6  8  7 ] ...
       [  1 11  2  5  4  3  8 ] ...
       [ 15 13 14  6  7 10  9 ] ];

YR = [ 'aaaaaab' 'bbbbbbb' 'cbbabaa' ...
       'ccccccc' 'cdcdccc' 'ddddddd' ...
       'dededdd' 'edeeeee' 'eeeeeee' ] - 'a' + 1;

% "Left" detector fibers

XL = [ [  7  6  4  5  8  7  5 ] ...
       [ 13  9 11 10 12  4  6 ] ...
       [ 10 11  7  8  9  6  8 ] ...
       [ 12  9 12 11 10  4  5 ] ...
       [  8  5  7  6  4 10 12 ] ...
       [  4  7  8  6  5 13  9 ] ...
       [ 12 13 11 10  9  7 11 ] ...
       [ 10 12 11 13  9  5  6 ] ...
       [  6  4  8  5  7  8  4 ] ];

YL = [ 'fffffgg' 'fffffgg' 'hhhhhhg' ...
       'iihiihh' 'jjjjjjj' 'iiiiijj' ...
       'gggggkj' 'kkkkkkk' 'lllllkk' ] - 'a' + 1;

X = [ XR(:); XL(:) ];
Y = [ YR(:); YL(:) ];

% Coordinates of user mouse-clicks

if (doRight)
   x = [  4 12  4 12 ]';
   y = [  1  1  5  5 ]';
else
   x = [  4 13  4 13 ]';
   y = [  6  6 11 11 ]';
end

% Compute mapping

[A,B] = calcAffine(x,y,u,v);    

% For each fiber, get its coordinates and build the mask

mask = cell(1,length(X));

for k = 1:length(X); 
   R = A * [ X(k) Y(k) ]' + B; 
   U(k) = R(1); 
   V(k) = R(2); 

   % Use flood-fill algorithm to find points at or above threshold
   [m,msk] = defiber(rawdata,U(k),V(k),threshold);

   % Save mask
   mask{k} = sparse(msk);
end

return;
