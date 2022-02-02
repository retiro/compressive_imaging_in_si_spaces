function opts = chkOptargs(opts)
%%CHKOPTARGS Checks validity of optional arguments.
%
% Copyright (c) 2013.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich.

% Set default values.
if ~isfield(opts, {'tol'}),      opts.tol       = 1.0e-8; end
if ~isfield(opts, {'tolpcg'}),   opts.tolpcg    = 1.0e-4; end
if ~isfield(opts, {'maxiters'}), opts.maxiters  = 100;    end
if ~isfield(opts, {'mxiterpcg'}),opts.mxiterpcg = 400;    end
if ~isfield(opts, {'maxstepsz'}),opts.maxstepsz = 0.999;  end
if ~isfield(opts, {'cnt1'}),     opts.cnt1      = 0.1;    end
if ~isfield(opts, {'cnt2'}),     opts.cnt2      = 0.5;    end
if ~isfield(opts, {'cnt3'}),     opts.cnt3      = 0.8;    end
if ~isfield(opts, {'ap'}),       opts.ap        = 0.5;    end
if ~isfield(opts, {'ac'}),       opts.ac        = 0.1;    end
if ~isfield(opts, {'linsol'}),   opts.linsol    = 0;      end
if ~isfield(opts, {'verbose'}),  opts.verbose   = 1;      end

% Check optional arguments.
if opts.verbose < 0, 
    error('The argument verbose must be non-negative.'); 
end
if opts.linsol ~= 0 && opts.linsol ~= 1, 
    error('The argument linsol must 0 or 1.'); 
end
if opts.maxstepsz <= 0 || opts.maxstepsz > 1, 
    error(['The argument maxstepsz must be greater than 0 and smaller or ',...
        'equal to 1.']); 
end
if opts.tol <= 0, 
    error('The argument tol must be greater than 0.'); 
end
if opts.tolpcg <= 0, 
    error('The argument tol must be greater than 0.'); 
end
if opts.maxiters <= 0, 
    error('The argument maxiters must be greater than 0.'); 
end
if opts.mxiterpcg <= 0, 
    error('The argument maxiters must be greater than 0.'); 
end
if opts.cnt1 <= 0 || opts.cnt1 > 1, 
    error(['The argument cnt1 must be greater than 0 and smaller or ',...
        'equal to 1.']); 
end
if opts.cnt2 <= 0 || opts.cnt2 > 1, 
    error(['The argument cnt2 must be greater than 0 and smaller or ',...
        'equal to 1.']); 
end
if opts.cnt3 <= 0 || opts.cnt3 > 1, 
    error(['The argument cnt3 must be greater than 0 and smaller or ',...
        'equal to 1.']); 
end
if opts.ap <= 0 || opts.ap >= 1,
error(['The argument ap must be greater than 0 and smaller ',...
       'than 1.']);
end
if opts.ac <= 0 || opts.ac >= 1,
error(['The argument ac must be greater than 0 and smaller ',...
       'than 1.']);
end
if opts.ac >= opts.ap,
error('The argument ac must be smaller than ap ');
end

end
