% MOVESRCSLAB  Move optodes 1 scattering length into the medium (slab)
%
% Src = moveSrcSlab(Src, Thickness, Muspo);
%
%   Src       Optode position vector
%   Thickness Thickness of slab (inf for semi-infinite slabs)
%   Muspo     Background scattering coefficient
%
% Returns:
%   Src       Displaced source vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pSrc] = moveSrcSlab(Src, Thickness, Muspo);

if (size(Src,2) ~= 3)
   error('Malformed source vector');
end

if (length(Muspo) ~= 1)
   error('Malformed mu_s')
end
   
for k = 1:size(Src,1)
   if     (Thickness == 0)
      error('Illegal zero thickness');
   elseif (Thickness > 0)
      % Move optodes from air to interface
      
      if (Src(k,3) < 0)
	 disp(sprintf('Warning: moving source %d from air to iterface',k));
	 disp(sprintf('\t[%f,%f,%f] -> [%f,%f,0.0]\n', ...
		      Src(k,1), Src(k,2), Src(k,3), Src(k,1), Src(k,2)));
	 
	 Src(k,3) = 0; 
      end
      
      if (~isinf(Thickness))
	 if (Src(k,3) > Thickness); 
	    disp(sprintf('Warning: moving source %d from air to iterface',k));
	    disp(sprintf('\t[%f,%f,%f] -> [%f,%f,%f]\n', ...
			 Src(k,1), Src(k,2), Src(k,3), ...
			 Src(k,1), Src(k,2), Thickness));
	    
	    Src(k,3) = Thickness; 
	 end
      else
	 if (Src(k,3) <         0); 
	    disp(sprintf('Warning: moving source %d from air to iterface',k));
	    disp(sprintf('\t[%f,%f,%f] -> [%f,%f,0.0]\n', ...
			 Src(k,1), Src(k,2), Src(k,3), Src(k,1), Src(k,2)));
	    
	    Src(k,3) =         0; 
	 end
      end
      
      % Move one scattering length into medium
      
      if (Src(k,3) < Thickness/2)
	 Src(k,3) = Src(k,3) + 1/Muspo;
      else
	 Src(k,3) = Src(k,3) - 1/Muspo;
      end
   elseif (Thickness < 0)
      % Move optodes from air to interface
      
      if (Src(:,3) > 0)
	 disp(sprintf('Warning: moving source %d from air to iterface',k));
	 disp(sprintf('\t[%f,%f,%f] -> [%f,%f,0.0]\n', ...
		      Src(k,1), Src(k,2), Src(k,3), Src(k,1), Src(k,2)));

	 Src(:,3) = 0; 
      end

      if (~isinf(Thickness))
	 if (Src(:,3) < Thickness)
	    disp(sprintf('Warning: moving source %d from air to iterface',k));
	    disp(sprintf('\t[%f,%f,%f] -> [%f,%f,%f]\n', ...
			 Src(k,1), Src(k,2), Src(k,3), ...
			 Src(k,1), Src(k,2), Thickness));

	    Src(:,3) = Thickness; 
	 end
      else
	 if (Src(:,3) >         0); 
	    disp(sprintf('Warning: moving source %d from air to iterface',k));
	    disp(sprintf('\t[%f,%f,%f] -> [%f,%f,0.0]\n', ...
			 Src(k,1), Src(k,2), Src(k,3), Src(k,1), Src(k,2)));

	    Src(:,3) =         0; 
	 end
      end
      
      % Move one length into medium
      
      if (Src(k,3) < Thickness/2)
	 Src(k,3) = Src(k,3) + 1/Muspo;
      else
	 Src(k,3) = Src(k,3) - 1/Muspo;
      end
   end
end

pSrc = Src;

return;
