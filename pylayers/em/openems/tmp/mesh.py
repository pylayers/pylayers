mport pdb
from xml.etree.ElementTree import Element
from xml.etree import ElementTree

%function [lines quality] = AutoSmoothMeshLines( lines, max_res, ratio, varargin)
% function [lines quality] = AutoSmoothMeshLines( lines, max_res, ratio, varargin)
%
%   Generate smooth mesh lines by choosing an appropriate algorithm.
%
%   Currently supported algorithm:
%       SmoothMeshLines, SmoothMeshLines2 and RecursiveSmoothMesh
%
%  arguments:
%   lines:      given fixed lines to create a smooth mesh in between
%   max_res:    desired max. resolution
%   ratio:      grading ratio: desired neighboring line-delta ratio
%                   - default is 1.5
%                   - see also 'allowed_max_ratio' argument
%
%  variable arguments ('keyword',value):
%   algorithm:          define subset of tried algorihm, e.g. [1 3]
%   symmetric:          0/1 force symmetric mesh (default is input symmetry)
%   homogeneous:        0/1 force homogeneous mesh
%   allowed_min_res:    allow a given min resolution only
%   allowed_max_ratio:  allow only a given max. grading ratio
%                           (default --> ratio*1.25)
%   debug:              0/1 off/on
%
% example:
%   lines = AutoSmoothMeshLines([-100 -10 10 100], 20, 1.5, 'algorihm', ...
%   1:3);
%
% See also InitCSX, DefineRectGrid
%
% CSXCAD matlab interface
% -----------------------
% author: Thorsten Liebig (C) 2012

if (nargin<2)
    error('CSXCAD:AutoSmoothMeshLines','lines and max_res are a required parameter');
end

if (nargin<3)
    ratio = 1.5;
end

lines = sort(unique(lines));

range = lines(end)-lines(1);
if (~isempty(find(diff(lines)<range*1e-6)))
    warning('CSXCAD:AutoSmoothMeshLines','some lines found with very small distance which may cause smoothing failure!');
end

methods = {};
methods{end+1} = @SmoothMeshLines;
methods{end+1} = @SmoothMeshLines2;
methods{end+1} = @RecursiveSmoothMesh;

requires_homogen = 0;
requires_symmetric = CheckSymmtricLines(lines);
allowed_min_res = 0;
debug = 0;
max_ratio = ratio*1.25;

algorithm = 1:numel(methods);

for vn=1:2:numel(varargin)
    if (strcmpi(varargin{vn},'algorithm'))
        algorithm = intersect(varargin{vn+1},algorithm);
    end
    if (strcmpi(varargin{vn},'symmetric'))
        requires_symmetric = varargin{vn+1};
    end
    if (strcmpi(varargin{vn},'homogeneous'))
        requires_homogen = varargin{vn+1};
    end
    if (strcmpi(varargin{vn},'force_min_res'))
        requires_homogen = varargin{vn+1};
    end
    if (strcmpi(varargin{vn},'allowed_min_res'))
        allowed_min_res = varargin{vn+1};
    end
    if (strcmpi(varargin{vn},'allowed_max_ratio'))
        max_ratio = varargin{vn+1};
    end
    if (strcmpi(varargin{vn},'debug'))
        debug = varargin{vn+1};
    end
end

for m=algorithm
    if (debug>0)
      disp(['AutoSmoothMeshLines: trying method: ' func2str(methods{m})]);
      tic
    end
    out_lines{m} = methods{m}(lines, max_res, ratio, 'CheckMesh', false);
    if (debug>0)
      toc
    end
    quality(m) = eval_mesh(out_lines{m}, max_res, max_ratio, requires_homogen, requires_symmetric, allowed_min_res, 1);
    if (quality(m)==100)           % uncomment to release!
       lines = out_lines{m};
       if (debug>0)
           disp(['AutoSmoothMeshLines: The winner with 100% is ' func2str(methods{m})]); % remove to release
       end
       return
    end
end
winner = find(quality==max(quality),1);
lines = out_lines{winner};
if (debug>0)
    disp(['New_SmoothMeshLines: The winner with ' num2str(quality(winner)) '% is ' func2str(methods{winner})]);  % remove to release
end

% show mesh problems
eval_mesh(lines, max_res, ratio, requires_homogen, requires_symmetric, allowed_min_res, 0);
return

error('CSXCAD:AutoSmoothMeshLines',['unknown algorithm requested: ' algorithm ]);

end


def eval_mesh(**kwargs):
    """
    """
    defaults = {'lines':,
                'max_res':,
                'ratio':,
                'requires_homogen':,
                'requires_symmetric':,
                'allowed_min_res':,
                'silent':}

    quality = 100

    results = AnalyseMesh(lines);
    if ((requires_homogen==1) && (results.homogeneous~=1))
    if (silent==0)
        warning('CSXCAD:AutoSmoothMeshLines','method failed to create homogenous mesh');
    end
    quality = -1;
    return
end
if ((requires_symmetric==1) && (results.symmetric~=1))
    if (silent==0)
        warning('CSXCAD:AutoSmoothMeshLines','method failed to create symmetric mesh');
    end
    quality = -1;
    return
end

if ((allowed_min_res>0) && (results.min_res<allowed_min_res))
    if (silent==0)
        warning('CSXCAD:AutoSmoothMeshLines','method failed to obey allowed min res!');
    end
    quality = -1;
    return
end

if (results.max_res>max_res*1.01)
    if (silent==0)
        warning('CSXCAD:AutoSmoothMeshLines',['method failed to fulfill max. res: ' num2str(results.max_res) ' > ' num2str(max_res)]);
    end
    quality = quality*(max_res/results.max_res);
end
if (results.max_ratio>ratio*1.01)
    if (silent==0)
        warning('CSXCAD:AutoSmoothMeshLines',['method failed to fulfill the max. ratio: ' num2str(results.max_ratio) ' > ' num2str(ratio)]');
    end
    quality = quality*(ratio/results.max_ratio);
end
end
