function testDocExamples()
%TESTDOCEXAMPLES  Every "Example:" block in the docstrings must run.
%   Extracts the indented code under each "Example" heading of every
%   matlab/src/*.m docstring and evaluates it in an isolated workspace.
%   Lines starting with "Note" or "See also" end the block.  The Python
%   port runs its docstring examples through pytest --doctest-modules.
end

%!function [ok, msg] = run_doc_example(src__)
%!  ok = true; msg = '';
%!  try
%!    evalc(src__);
%!  catch err__
%!    ok = false; msg = err__.message;
%!  end

%!test
%! src = fileparts(which('elliptic12'));   % mfilename is empty inside test blocks under test(<full path>)
%! files = dir(fullfile(src, '*.m'));
%! nrun = 0;
%! for f = 1:numel(files)
%!   lines = strsplit(fileread(fullfile(files(f).folder, files(f).name)), "\n");
%!   i = 1;
%!   while i <= numel(lines)
%!     if ~isempty(regexp(lines{i}, '^\s*%\s*Example', 'once'))
%!       code = {}; j = i + 1;
%!       while j <= numel(lines) && ~isempty(regexp(lines{j}, '^\s*%\s{3,}\S', 'once'))
%!         c = regexprep(lines{j}, '^\s*%\s*', '');
%!         if ~isempty(regexp(c, '^\s*(Note|See also)', 'once')), break; end
%!         if isempty(regexp(c, '^\s*%', 'once')), code{end+1} = c; end
%!         j = j + 1;
%!       end
%!       if ~isempty(code)
%!         nrun = nrun + 1;
%!         [ok, msg] = run_doc_example(strjoin(code, "\n"));
%!         assert(ok, sprintf('%s:%d docstring example failed: %s', files(f).name, i, msg));
%!       end
%!       i = j;
%!     else
%!       i = i + 1;
%!     end
%!   end
%! end
%! assert(nrun >= 8, sprintf('expected at least 8 docstring examples, found %d', nrun));

