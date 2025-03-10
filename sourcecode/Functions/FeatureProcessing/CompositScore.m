function Score = CompositScore(x,y)
arguments 
    x (:,1) {mustBeNumeric,mustBeReal}
    y (:,1) {mustBeNumeric,mustBeReal,mustBeEqualSize(x,y)}
end
C = 999;

Nu = sum(x~=0); %number of non zero elements in x
Nlu = sum(x~=0 & y~=0); %number of non zero elements in both vectors

% calculate cosine similarity
matchFactor = dot(x,y)^2/(sum(x.^2)*sum(y.^2));

%remove zero elements
idx = x==0 | y==0;

x(idx) = [];
y(idx) = [];

% shift vectors by 1
downX = [0;x];
downY = [0;y];
x = [x;0];
y = [y;0];

% calculate composit
composit = (x./downX.*downY./y);
composit(~isfinite(composit)) = [];
composit(composit>1) = composit(composit>1).^-1;

Fr = 1/Nlu * sum(composit);

%% calculate composite score
Score = C*(Nu*matchFactor + Nlu*Fr)/(Nu + Nlu);
Score = round(Score);
end

% Custom validation function
function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end