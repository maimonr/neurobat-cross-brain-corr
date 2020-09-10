function resultantSet = multiIntersect(varargin)

if nargin < 2
    resultantSet = varargin{1};
    return
end

resultantSet = intersect(varargin{1},varargin{2});

for set_k = 3:nargin
    resultantSet = intersect(resultantSet,varargin{set_k});
end

end