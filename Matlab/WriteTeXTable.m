function [] = WriteTeXTable(fid, header, style, tableData, above_tabular, below_tabular)
%% WriteTeXTable - Writing tex tables from Matlab.
%
% NOTE: below, I will use NCOL and NROW to denote the maximal number of
% columns in the TeX table ('maximal' in the sense that I ignore
% multicolumn), while NROW is the number of rows of data, excluding
% header information.
%
% Input arguments
% ---------------
% FID         File ID you got before running this program by executing
%             something like
%
%               fid = fopen('TableName.tex', 'w')
%
% HEADER      Either:
%             1.  Cell array of column titles/headers
%             2.  An array structure of length Nheadrows, each with
%                 field 'header' that contains a cell array of column
%                 titles/headers
%
%             NOTE: A cell array of column titles/headers can look like
%             this
%
%                  {'Col1', 'Col2', 'Col3', 'Col4'};
%             OR
%                  {'Multi1', [1 2], 'Col3', 'Col4'};
%
%             where having an array after text will cause that text to
%             span the columns specified. Therefore, only text entries
%             show up as column titles/headers, while arrays can be used
%             to signify multi-column spanning
%
% STYLE       Latex style for the columns like 'r|cccc' or 'l|rr|rr|'
%
% TABLEDATA   Array or cell matrix of table content to write. Though
%             headers are treated differently from data to display, no
%             such distinction is made for the first column or couple of
%             columns for labels. Just throw them in here.
%
%             If TABLEDATA is a numeric array, the numbers will be
%             displayed to significant digits determined by this
%
%             Text is treated as given (i.e. it should be written as
%             is). Numeric entries in 
%
%             Example table_data
%
%               table_data = [{'row1'; 'row2'; 'row3'},...
%                             num2cell(magic(3))];
%
% caption     Optional caption for table
% fontsize    Optional fontsize like \scriptsize or \footnotesize
%
% Example usage
% -------------
%
%   row_names  = {'row1'; 'row2'; 'row3'};
%   data2write = randn(3);
%   header     = {'col1', 'col2', 'col3'};
%   table_data = [row_names, num2cell(data2write)];
%   style      = 'r|ccc';
%
%   fid = fopen('file.tex', 'w');
%   write_table_tex(fid, header, style, table_data, 'Random Numbers')
%   fclose(fid);
%

Nrow = size(table_data,1);
Ncol = size(table_data,2);

% Can set fileid to NaN to print to screen
if isnan(fid)
  writeline = @(s) fprintf('%s\n', s); % Might be defined in calling function
else
  writeline = @(s) fprintf(fid, '%s\n', s); % Might be defined in calling function
end
%helper_functions; 
iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

%% Write the tex part with the table sizing
writeline('\begin{table}[htpb!]');
if exist('fontsize', 'var')
  writeline(fontsize);
end
writeline('\centering');
writeline(['\begin{tabular}{', style, '}']);


%% Write header 

  % Handle an extra grouping in the header 
  if sum(cellfun(@(c) iscell(c), header)) 
    % ^If the "header" cell has cell elements, then the header will be multi-line

    % How many columns the multicolumn should be
    col_width = cellfun(@(h) iif(iscell(h), @() length(h{2}), true, 1), header)';
      % ^Need the extra @() length(h{2}) in there so we don't evaluate h{2} right
      % away, which would throw an error if h is not a cell

    % Construct the multicolumn tex headers
    ismulti = (col_width > 1);
    multi_headers          = repmat({''}, 1, length(ismulti));
    multi_headers(ismulti) = ...
      arrayfun(@(h) sprintf('\\multicolumn{%d}{c}{%s}', col_width(h), header{h}{1}), ...
               find(ismulti), 'un', 0);

    % Write the multicolumn tex headers
    write_line = writeline([strjoin(multi_headers, ' & '), '\\']);

    % Write the clines
    clines = [cumsum(col_width)-col_width+1, cumsum(col_width)];
    clines = clines(ismulti,:);
    arrayfun(@(ln) writeline(sprintf('\\cline{%d-%d}', clines(ln,1), clines(ln,2))), 1:size(clines,1));

    % Make the header into the last and final header row (the one closest to the data)
    header = cellfun(@(h) iif(iscell(h), @() h{2}, true, h), header, 'un', 0);
    header = vertcat(header{:})';
  end

  % Write the usual header; if there is a column with NaN, extend the
  % column to the left to cover it via multicolumn. i.e. 
  %
  % {'Jerry', 'George', NaN, 'Kramer', NaN, NaN, 'Elaine'} %
  % will be translated as
  %
  % Jerry & \multicolumn{2}{c}{George} & \multicolumn{3}{c}{Kramer} & Elaine\\
  colspans = cumsum(~cellfun(@(h) isnumeric(h) && isnan(h), header)); 
  nlabs = length(unique(colspans)); % number of header labels that aren't empty
  fullheader = cell(nlabs,1);
  for c = 1:nlabs
    span = sum(colspans == c);
    lab  = header{find(colspans == c, 1)};
    if span-1
      fullheader{c} = sprintf('\\multicolumn{%d}{c}{%s}', span, lab);
    else
      fullheader{c} = lab; % Account for extra slash as escape
    end
  end
  if size(fullheader,1) > 1 && size(fullheader,2) == 1, fullheader = fullheader'; end
  writeline([strjoin(fullheader, ' & '), '\\\hline\hline']);


%% Write the table

  % Function to format the table entries/elements
  fmt = @(x) iif(isnumeric(x) && abs(x) < 1,  sprintf('%.3f', x), ...
                  isnumeric(x) && abs(x) < 10,  sprintf('%.2f', x), ...
                  isnumeric(x),                 sprintf('%.1f', x), ...
                  true,                         x);

  % Loop over rows, possibly extending using multicolumn like above
  for r = 1:size(table_data,1)

    row = table_data(r,:); 

    % Array that will mark different column spans; will look like 
    %   [1 2 2 3 4 4 4] 
    % denoting that the first entry in "fullrow" below will span 1 column, the
    % second will span 2, the third will span 1, and the fourth will span 3
    % columns

    colspans = cumsum(~cellfun(@(h) isnumeric(h) && isnan(h), row)); 
    nentries = length(unique(colspans)); % number of entries (accounting for any multicolumn business)

    fullrow = repmat({''}, nentries,1); % Will be a cell with the number of unique columns entries to be filled 
    for c = 1:nentries
      span  = sum(colspans == c); % Number of columns the first entry must span
      entry = row{find(colspans == c, 1)}; % Get that first entry
      if span-1
        fullrow{c} = sprintf('\\multicolumn{%d}{c}{%s}', span, fmt(entry));
      else
        fullrow{c} = fmt(entry);
      end
    end
    if size(fullrow,1) > 1 && size(fullrow,2) == 1, fullrow = fullrow'; end
    writeline([strjoin(fullrow, ' & '), '\\']);
  end


% Conclude
writeline('\hline');
writeline('\end{tabular}');

% Caption
if exist('caption', 'var')
  writeline(sprintf('\\caption{%s}', caption));
end
writeline('\end{table}');

writeline('');

end
