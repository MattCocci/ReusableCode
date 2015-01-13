%% OVERVIEW - write_figs_tex.m
%
% This is a recursive function to write \includefigure statements for a bunch
% of plots.
%
% This function takes a nested cell with information about where plots are
% saved, what their captions are, and how the plots should be structured
% (sections, subsections, subsubsections). It then writes all the tex code
% necessary to include those figures in a document, preserving the structure,
% writing sections and section names, and writing captions.
%
% STRUCTURE OF CELL
% -----------------
% Since this function is recursive, the cell should be set up with nested
% versions of the same pattern.
%  
% Say that we want to create a document that has multiple pictures of
% characters in various sitcoms.  The document might look like this:
% 
% \section{Seinfield}
% 
% \subsection{Jerry}
% figure1
% figure2
% 
% \subsection{Elaine}
% figure1
% figure2
% figure3
% 
% \subsection{George}
% figure1
% figure2
% 
% \subsection{Kramer}
% figure1
%
% 


function [ ] = write_figcell_tex(fid, obj, n)

  if n > 3
    error('Too many nestings in the cell object. Only three nestings allowed since LaTeX only defined up to subsubsection.')
  end

  if ~exist('writeline')
    writeline = @(s) fprintf(fid, '%s\n', s); % Might be defined in calling function
  end


  Nsections = size(obj, 1);
  for s = 1:Nsections
    row = obj(s,:);

    % Write section header
    headers = {'\section', '\subsection', '\subsubsection'};
    writeline(sprintf('\n\n%s\n%s{%s}', '\clearpage', headers{n}, row{1}));

    % If all strings, write them
    if sum(sum(cellfun(@(s) ~isstr(s), row{2}))) == 0
      captions = row{2}(:,1);
      files    = row{2}(:,2);
      Nfiles   = length(files);
      for f = 1:Nfiles
        writeline('');
        writeline('\begin{figure}[h!]');
        writeline('\centering');
        writeline(sprintf('\t%s{%s}', '\caption', captions{f}));
        writeline(sprintf('\t%s{%s}', '\includegraphics[scale=0.6,trim={2.5cm, 6, 1.7cm, 7cm}, clip]', files{f}));
        writeline('\end{figure}');
      end

    else
      write_figcell_tex(row{2}, n+1);
    end

  end


end
