function uniformAxisSpacing(ax,margins,spacing)

margins = fillEmptyMargins(margins,ax);
spacing = fillEmptySpacing(spacing,ax);

ncols = size(ax,2);
xstarts = zeros(size(ax));
xstarts(:,1) = margins.left;

width_panel_region = 1 - margins.left - margins.right;
axis_width = (width_panel_region - (ncols-1)*spacing.horizontal)/ncols;
for ci = 2:ncols
    xstarts(:,ci) = xstarts(:,ci-1) + axis_width + spacing.horizontal;
end

nrows = size(ax,1);
ystarts = zeros(size(ax));
ystarts(end,:) = margins.bottom;

height_panel_region = 1 - margins.top - margins.bottom;
axis_height = (height_panel_region - (nrows-1)*spacing.vertical)/nrows;
for ri = (nrows-1):-1:1
    ystarts(ri,:) = ystarts(ri+1,:) + axis_height + spacing.vertical;
end

for ri = 1:nrows
    for ci = 1:ncols
        ax(ri,ci).Position = [xstarts(ri,ci),ystarts(ri,ci),axis_width,axis_height];
    end
end

end

function margins = fillEmptyMargins(margins,ax)

if ~isfield(margins,"left") || isempty(margins.left)
    margins.left = ax(1,1).Position(1);
end

if ~isfield(margins,"right") || isempty(margins.right)
    margins.right = 1-ax(1,end).Position(1)-ax(1,end).Position(3);
end

if ~isfield(margins,"bottom") || isempty(margins.bottom)
    margins.bottom = ax(end,1).Position(2);
end

if ~isfield(margins,"top") || isempty(margins.top)
    margins.top = 1-ax(1,1).Position(2)-ax(1,1).Position(4);
end

end

function spacing = fillEmptySpacing(spacing,ax)

if numel(ax)==1
    spacing.horizontal = 0;
    spacing.vertical = 0;
    return
end

if size(ax,1)==1
    spacing.vertical=0;
elseif ~isfield(spacing,"vertical") || isempty(spacing.vertical)
    spacing.vertical = ax(1,1).Position(2) - (ax(2,1).Position(2)+ax(2,1).Position(4));
end


if size(ax,2)==1
    spacing.horizontal=0;
elseif ~isfield(spacing,"horizontal") || isempty(spacing.horizontal)
    spacing.horizontal = ax(1,2).Position(1) - (ax(1,1).Position(1)+ax(1,1).Position(3));
end

end
