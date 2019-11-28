function dg = torGradient(func, points, step_size, varargin)
    center = func(points, varargin{1}, varargin{2});
    dg = zeros(size(points));
    for i=1:size(points,1)
        dir_multiplier = zeros(size(points,1),1);
        dir_multiplier(i) = 1;
        dg(i,:) = ((func(points + (dir_multiplier .* step_size), varargin{1}, varargin{2}) - center) + (center - func(points - (dir_multiplier .* step_size), varargin{1}, varargin{2})))/(2*step_size);
    end
end