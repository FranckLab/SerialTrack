function I = seedBeadsN(sigma, OS, sizeI)


[nBeads, nDims] = size(OS);
if nargin < 3, sizeI = ceil(max(OS) - min(OS) + 1); end
if numel(sigma) == 1, sigma = sigma*ones(1,nDims); end

beadSize = ceil(sigma)*2*5;

if nDims == 2,
    [m{1}, m{2}] = ndgrid(1:beadSize(1),1:beadSize(2));
else
    [m{1}, m{2}, m{3}] = ndgrid(1:beadSize(1),1:beadSize(2),1:beadSize(3));
end

OS_ = bsxfun(@plus, OS - floor(OS), beadSize/2);

%%
f = num2cell(zeros(nBeads,1));
for i = 1:nBeads
    for j = 1:nDims
        f{i} = f{i} - ( (m{j} - OS_(i,j)) / (2*sigma(j))).^2;
    end
    f{i} = exp(f{i});
end

idx = cell(nDims,1);
mask = cell(nDims,1);

for i = 1:nDims
    idx{i} = bsxfun(@plus, floor(OS(:,i)), repmat(1:beadSize(i),[nBeads,1])) - beadSize(i)/2;
    mask{i} = idx{i} < 1 | idx{i} > sizeI(i);
    idx{i} = mat2cell(idx{i}, ones(1,nBeads),beadSize(i));
    
    j_ = find(any(mask{i},2));
    if ~isempty(j_)
        for j = j_(1):j_(end)
            switch i
                case 1, f{j}(mask{1}(j,:),:,:) = [];
                case 2, f{j}(:,mask{2}(j,:),:) = [];
                case 3, f{j}(:,:,mask{3}(j,:)) = [];
            end
            
            idx{i}{j}(mask{i}(j,:)) = [];
        end
    end
    
end


I = zeros(sizeI);
if nDims ==2
    
    for i = 1:nBeads
        I(idx{1}{i},idx{2}{i}) = I(idx{1}{i},idx{2}{i}) + f{i};
    end
    
else
    
    for i = 1:nBeads
        I(idx{1}{i},idx{2}{i},idx{3}{i}) = I(idx{1}{i},idx{2}{i},idx{3}{i}) + f{i};
    end
    
end

I(I > 1) = 1;


end