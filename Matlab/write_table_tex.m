function [] = write_table_tex(f, header, style, table_data)

Ncol  = length(table_data);
Nrow  = cellfun(@(c) length(c), table_data);
if ~all(Nrow == Nrow(1))
  error('The number of rows in the table differs across columns')
else
  Nrow = Nrow(1);
end

writeline = @(s) fprintf(f, '%s\n', s); % Might be defined in calling function
helper_functions;


%% Write the tex part with the table sizing
writeline('\begin{table}')
writeline('\centering')
writeline(['\begin{tabular}{', style, '}'])


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

  % Write the usual header
  writeline([strjoin(header, ' & '), '\\\hline\hline']);


%% Write the table

  % Function to format the table entries/elements
  fmt = @(x) iif(isnumeric(x) && abs(x) < 1,  sprintf('%9.3f', x), ...
                isnumeric(x) && abs(x) < 10,  sprintf('%9.2f', x), ...
                isnumeric(x),                 sprintf('%9.1f', x), ...
                true,                         x);


  % Function to write a row
  getrow   = @(r) arrayfun(@(c) fmt(table_data{c}{r}), 1:Ncol, 'un', 0);
  writerow = @(r) writeline([strjoin(getrow(r), ' & '), '\\']);

  % Write all the rows
  arrayfun(@(r) writerow(r), 1:Nrow);


% Conclude
writeline('\end{tabular}');
writeline('\end{table}');
writeline('');


end
