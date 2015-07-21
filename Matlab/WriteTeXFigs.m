function [] = WriteTeXFigs(fid, captions, files, sizing)
%% writeFigsTex - For writing \includegraphics to tex docs

  % Make sure files and captions are stored in cells for looping below
  if ischar(files)
    files = {files};
  end
  if ischar(captions)
    captions = {captions}; % Oh caption, my caption
  end
  Nfiles = length(files);

  % Set figure sizing param; if not given, provide default; if there is
  % one, but it's less than the number of files, use it for all files
  if ~exist('sizing', 'var')
    sizing = 'scale=0.5,trim={2.0cm, 5cm, 1.7cm, 5.75cm}, clip';
  end
  if ischar(sizing)
    sizing = {sizing};
  end
  if length(sizing) == 1 && Nfiles > 1
    sizing = repmat(sizing, Nfiles, 1);
  end


  % Function to write things can write writeline('\includegraphics')
  % this way, rather than fprintf('\\includegraphics') or something
  if isempty(fid) || isnan(fid)
    writeline = @(s) fprintf('%s\n', s); % For debugging, print to screen
  else
    writeline = @(s) fprintf(fid, '%s\n', s);
  end


  % Loop and write
  for n = 1:Nfiles
    if n > 1 && mod(n,2) == 1
      writeline('\clearpage');
    end
    writeline('\begin{figure}[htpb!]');
    writeline('\centering');
    writeline(sprintf('  %s{%s}', '\caption', captions{n}));
    writeline(sprintf('  %s[%s]{%s}', '\includegraphics', sizing{n}, files{n}));
    writeline('\end{figure}');
  end

end
