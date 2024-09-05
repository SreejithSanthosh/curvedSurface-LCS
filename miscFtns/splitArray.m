function splitCell = splitArray(a,n)
a = reshape(a,[numel(a),1]);
b = mod(numel(a),n);
X = fix(numel(a)/n);

if b == 0
    Y = X*ones(1,n);
else 
    Y = [X*ones(1,n-1),X+b];
end 

% Y = [n*ones(1,fix(X/n)),1+mod(X,n)];

splitCell = mat2cell(a,Y,1);
end 