function f = figureOnRight(varargin)

f = figure(varargin{:});
f.Units = "normalized";
f.Position(1) = 1 - f.Position(3);