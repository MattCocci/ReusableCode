function [] = write_preamble(fid, title, extra)

  writeline = @(str) fprintf(fid, '%s\n', str);

  % Set up the preamble
  writeline('\documentclass[11pt]{article}');
  writeline(sprintf('%s{%s}', '\title', title));
  writeline('\author{Matthew Cocci}');
  writeline('\date{\today}');
  writeline('\usepackage{amsmath}');
  writeline('\usepackage{graphicx}');
  writeline('\usepackage{pdfpages}');
  writeline('');
  writeline('\usepackage[margin=1in]{geometry}');
  writeline('\usepackage{graphicx}');
  writeline('\usepackage{subfigure}');
  writeline('');

  writeline('\usepackage{hyperref}');
  writeline('\hypersetup{');
  writeline(' colorlinks,');
  writeline(' citecolor=black,');
  writeline(' filecolor=black,');
  writeline(' linkcolor=black,');
  writeline(' urlcolor=black ');
  writeline('}');

  cellfun(writeline, extra);

  writeline('\begin{document}');

end
