% LOADIMAGE  read in binary image from TD software
%
% [data,hdr] = loadimage(filename[,maxframe]);

function[data, hdr] = loadimage2(name, maxframe)

if (~exist('maxframe','var') | maxframe < 1)
   maxframe = inf;
end

% Open the file read-only

fid = fopen(name, 'rb');

if (fid < 0)
   error(['Error opening file ' name]);
end

% Figure out what version file this is

[mark,m] = fread(fid, 1, 'uint16');

if (m < 1)
   fclose(fid);
   error(['Error reading from file ' name ]);
else
   frewind(fid);
end

if     (mark == hex2dec('DDCC')) % low word
   versn = 2;
elseif (mark == hex2dec('FFFF'))
   versn = 1;
else
   % Obsolete
   versn = 0;
end

% Read the file header and the file data

wh = [];

if (versn == 2)
   % New-style file, stream of records
   iFrame = 1;
   
   while (~feof(fid) & (iFrame <= maxframe))
      tmphdr = getNewHdr(fid);
      x = tmphdr.x;
      y = tmphdr.y;
      
      if (isempty(x) & isempty(y))
         % Hit end of file, but not picked up by feof().
         % Break out of loop now
         
         close(wh);
         wh = [];

         break;
      else
         hdr(iFrame) = tmphdr;
         clear tmphdr;
      end
      
      [rawdata,ndata] = fread(fid, x*y, 'uint16');
      
      if (ndata ~= x * y)
         close(wh);
         wh = [];

         error(['Error reading from file ' name '; unexpected EOF']);
      end
      
      % Turn unsigned to signes shorts, only using first 12 bits anyway
      rawdata = reshape(int16(rawdata), x, y);
      
      data(:,:,iFrame) = rawdata;
      
      if (iFrame == 1)
         here = ftell(fid);
         fseek(fid, 0, 'eof');
         fend = ftell(fid);
         
         % Restore original file position
         fseek(fid, here, 'bof');
         
         % Guess how many frames of data are left
         nframes = floor(fend / here);
         
         % Do the smaller of the total file and the requested limit
         nframes = min(maxframe, nframes);
	 
         % Extend the data array
         if (nframes > 1)
           disp([ 'Estimating ' num2str(nframes) ' frames of data' ]);
           data(:, :, nframes) = 0;
         end
        
         clear here fend nframes;
      end
      
      if (size(data,3) > 1)
         if (isempty(wh))
            wh = waitbar(0.0, 'Loading data');
	    lastw = 0;
	    
            if (isempty(get(wh,'XDisplay')))
               % No X display, not a real box
               close(wh);
               wh = [];
            end
         end

         if (~isempty(wh))
	    % Update display in ~1% increments, not every frame
	    if (floor(100 * iFrame / size(data,3)) > lastw)
	       lastw = floor(100 * iFrame / size(data,3));
	       waitbar(lastw / 100, wh); drawnow;
	    end
         else
            disp(iFrame);
         end
      end

      iFrame = iFrame + 1;
   end

   if (~isempty(wh))
      close(wh);
   end
elseif (versn == 1)
   % Old-style file, 6-byte header, single frame
   
   hdr = getOldHdr(fid);
   
   [data, m] = fread(fid, 'uint16');     % Rest of file

   if (m ~= hdr.x * hdr.y)
      error(['Error reading from file ' name '; unexpected EOF']);
   end

   data = reshape(int16(data), hdr.x, hdr.y);
else
   fclose(fid);
   error('Unknown data header version');
end

% Done

fclose(fid);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[hdr] = getOldHdr(fid)

[foo, n] = fread(fid, 3, 'uint16');  % Header, if any

if (foo(1) == hex2dec('FFFF'))       % Got a valid header, set x,y
 	hdr.x = foo(2);
   hdr.y = foo(3);
else
   fclose(fid);
   error('getOldHdr called, no V1 header found');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[hdr] = getNewHdr(fid)

hdr.mark = fread(fid, 1, 'uint32');

if (hdr.mark ~= hex2dec('FFEEDDCC'))
   fclose(fid);
   error('getNewHdr() called, no V2 header found');
end

% Timestamp bytes are:
%  wYear, wMonth, wDayOfWeek, wDay, wHour, wMinute, wSecond, wMillisec
% I could repack it into a string, but I'm lazy today.
hdr.tstamp = fread(fid, 8, 'int16');

hdr.exposure = fread(fid, 1, 'uint32');
hdr.x        = fread(fid, 1, 'int16');
hdr.y        = fread(fid, 1, 'int16');

hdr.delay    = fread(fid, 1, 'int32');
hdr.lambda   = fread(fid, 1, 'int32');
hdr.muxx     = fread(fid, 1, 'int16');
hdr.muxy     = fread(fid, 1, 'int16');
hdr.imux     = fread(fid, 1, 'int32');

% Skip the reserved bytes
reserved     = fread(fid, 20, 'uint8');

return;
